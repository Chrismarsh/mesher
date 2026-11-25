"""Write the global mesh and parameter files read by CHM and its partition tool."""

import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import traceback

import numpy as np

if __package__ in (None, ''):  # launched by the configured MPI Python command
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mesher.mesh_store import batches, copy_array, disk_array, rank_face_range


def require_h5py():
    try:
        import h5py
        if not hasattr(h5py, 'File'):
            raise ImportError('h5py installation is missing the top-level Python package')
    except ImportError as exc:
        raise RuntimeError('write_hdf5 requires h5py. Install MPI-enabled h5py for parallel writes.') from exc
    return h5py


def _open_inputs(plan):
    geometry = {key: np.load(path, mmap_mode='r') for key, path in plan['geometry'].items()}
    params = {key: np.load(path, mmap_mode='r') for key, path in plan['parameters'].items()}
    return geometry, params


def _write_slices(plan, arrays, params, h5py, comm=None):
    rank, workers = (0, 1) if comm is None else (comm.rank, comm.size)
    driver = dict(driver='mpio', comm=comm) if workers > 1 else {}
    nfaces, nverts = len(arrays['elem']), len(arrays['vertex'])
    directory = Path(plan['directory'])

    # CHM reads rank-one datasets whose HDF5 datatype is a three-element array.
    # It reads the slash-prefixed attributes from the file root.
    with h5py.File(directory / 'mesh.h5', 'w', **driver) as mesh, \
            h5py.File(directory / 'param.h5', 'w', **driver) as parameters:
        # Parallel HDF5 requires all writers to create the same groups,
        # datasets and attributes in the same order, even on ranks whose
        # slice is empty. Rank-specific writes start after schema creation.
        group = mesh.create_group('mesh')
        vertices = group.create_dataset('vertex', (nverts,), dtype=np.dtype(('f8', (3,))))
        elements = group.create_dataset('elem', (nfaces,), dtype=np.dtype(('i4', (3,))))
        neighbours = group.create_dataset('neighbor', (nfaces,), dtype=elements.dtype)
        global_ids = group.create_dataset('cell_global_id', (nfaces,), dtype='>i4')
        sizes = group.create_dataset('local_sizes', (len(plan['sizes']),), dtype='>i4')
        permuted = len(sizes) > 1
        for name, value in [('proj4', plan['proj4']), ('version', '2.0.0' if permuted else '1.0.0'),
                            ('partition_method', 'metis' if permuted else '')]:
            mesh.attrs.create('/mesh/' + name, [value.encode('utf-8')], dtype='S256')
        mesh.attrs.create('/mesh/is_partition', [0], dtype='u1')
        mesh.attrs.create('/mesh/is_geographic', [plan['geographic']], dtype='u1')
        param_group = parameters.create_group('parameters')
        datasets = {key: param_group.create_dataset(key, (nfaces,), dtype='f8') for key in params}

        if rank == 0:
            sizes[:] = plan['sizes']
        lo, hi = rank_face_range(nverts, rank, workers)
        for start, end in batches(hi-lo):
            rows = slice(lo+start, lo+end)
            vertices[rows] = arrays['vertex'][rows]

        lo, hi = rank_face_range(nfaces, rank, workers)
        for start, end in batches(hi-lo):
            rows = slice(lo+start, lo+end)
            # CHM indexes parameter rows and neighbours by output face ID.
            # order finds the source face for each row. inverse translates
            # its neighbour IDs back to output rows, preserving boundary -1.
            source = arrays['order'][rows] if permuted else rows
            elements[rows] = arrays['elem'][source].astype('i4')
            adjacent = arrays['neigh'][source]
            if permuted:
                valid = adjacent >= 0
                adjacent[valid] = arrays['inverse'][adjacent[valid]]
            neighbours[rows] = adjacent.astype('i4')
            global_ids[rows] = np.arange(rows.start, rows.stop, dtype='i4')
            for key, dataset in datasets.items():
                dataset[rows] = params[key][source]


def export_hdf5(prefix, mesh, params, proj4, partitions, workers, launcher):
    nfaces, nverts = len(mesh['elem']), len(mesh['vertex'])
    if max(nfaces, nverts) > np.iinfo(np.int32).max:
        raise ValueError('CHM HDF5 requires face and vertex counts within signed 32-bit range')
    if not 1 <= partitions <= nfaces:
        raise ValueError(f'hdf5_partition_ranks must be between 1 and the number of faces ({nfaces})')
    if not proj4 or len(proj4.encode('utf-8')) >= 256:
        raise ValueError('CHM HDF5 requires a nonempty PROJ string under 256 bytes')
    for key, values in params.items():
        if not key or key == '.' or '/' in key or '\0' in key or values.shape != (nfaces,):
            raise ValueError(f'Invalid CHM HDF5 parameter {key!r}; expected one named column of {nfaces} faces')

    prefix = Path(prefix).absolute()
    mesh_path = prefix.with_name(prefix.name + '_mesh.h5')
    param_path = prefix.with_name(prefix.name + '_param.h5')
    with tempfile.TemporaryDirectory(prefix='.mesher-hdf5-', dir=prefix.parent) as work:
        directory = Path(work)
        # Give writers paths to the prepared arrays so each rank can map its
        # input on shared storage without receiving a copy of the full mesh.
        geometry = {key: os.fspath(mesh[key].filename) for key in ('vertex', 'elem', 'neigh')}
        sizes = [nfaces]
        if partitions > 1:
            from mesher.tools.meshpermutation import metis_face_order
            print(f'Ordering {nfaces} faces for {partitions} CHM ranks', flush=True)
            counts, order = metis_face_order(mesh['neigh'], partitions)
            sizes = counts.tolist()
            geometry.update(order=str(directory / 'order.npy'), inverse=str(directory / 'inverse.npy'))
            copy_array(geometry['order'], order, np.int32)
            inverse = disk_array(geometry['inverse'], (nfaces,), np.int32)
            for lo, hi in batches(nfaces):
                inverse[order[lo:hi]] = np.arange(lo, hi, dtype=np.int32)
            inverse.flush()
            del inverse, order

        plan = dict(geometry=geometry, parameters={k: os.fspath(v.filename) for k, v in params.items()},
                    sizes=sizes, proj4=proj4, geographic=bool(mesh['is_geographic']),
                    directory=work)
        print(f'Writing {mesh_path.name} and {param_path.name} with {workers} HDF5 writer(s)', flush=True)
        if workers == 1:
            _write_slices(plan, *_open_inputs(plan), require_h5py())
        else:
            manifest = directory / 'inputs.json'
            manifest.write_text(json.dumps(plan), encoding='utf-8')
            launch_writers(str(manifest), workers, launcher)

        # Publish the mesh last. An interrupted replacement must never pair a
        # previous mesh with parameters written in a new permutation.
        if mesh_path.is_dir():
            raise IsADirectoryError(mesh_path)
        if mesh_path.exists():
            os.replace(mesh_path, directory / 'previous_mesh.h5')
        try:
            os.replace(directory / 'param.h5', param_path)
        except OSError:
            if (directory / 'previous_mesh.h5').exists():
                os.replace(directory / 'previous_mesh.h5', mesh_path)
            raise
        os.replace(directory / 'mesh.h5', mesh_path)


def launch_writers(manifest, workers, launcher):
    from mpi4py import MPI

    args = [os.path.abspath(__file__), manifest, str(workers)]
    if launcher is not None:
        subprocess.check_call(f'{launcher} {shlex.join(args)}', shell=True)
    else:
        # Use shared storage for the result, as the other spawned Mesher stages
        # do. This also works when the launcher only has a control connection.
        shared = Path.cwd() if manifest == '--check' else Path(manifest).parent
        with tempfile.TemporaryDirectory(dir=shared, prefix='.hdf5-writers-') as work:
            result = Path(work) / 'result.json'
            child = MPI.COMM_SELF.Spawn(sys.executable, args=[*args, str(result)], maxprocs=workers)
            child.Disconnect()
            if not result.is_file():
                raise RuntimeError('HDF5 writers exited without completing the export')
            error = json.loads(result.read_text(encoding='utf-8'))
            if error:
                raise RuntimeError(error)


def _worker(manifest, expected, result=None):
    from mpi4py import MPI

    world = MPI.COMM_WORLD
    parent = MPI.Comm.Get_parent()
    # Agree on rank-local setup failures before anyone enters H5Fcreate.
    try:
        h5py = require_h5py()
        if world.size != expected:
            raise ValueError(f'HDF5 export expected {expected} writers, got {world.size}')
        if not h5py.get_config().mpi:
            raise RuntimeError('Parallel HDF5 export requires MPI-enabled h5py and HDF5; '
                               'set hdf5_nworkers=1 for a serial, bounded-memory write')
        if manifest != '--check':
            plan = json.loads(Path(manifest).read_text(encoding='utf-8'))
            arrays, params = _open_inputs(plan)
        error = None
    except Exception:
        error = traceback.format_exc()
    error = next((msg for msg in world.allgather(error) if msg), None)
    if not error and manifest != '--check':
        try:
            _write_slices(plan, arrays, params, h5py, world)
        except Exception:
            error = traceback.format_exc()
        error = next((msg for msg in world.allgather(error) if msg), None)
    if parent != MPI.COMM_NULL:
        try:
            if world.rank == 0:
                Path(result).write_text(json.dumps(error), encoding='utf-8')
        finally:
            parent.Disconnect()
        return
    if error:
        if world.rank == 0:
            print(error, file=sys.stderr, flush=True)
        sys.exit(1)


if __name__ == '__main__':
    _worker(sys.argv[1], int(sys.argv[2]), *sys.argv[3:])
