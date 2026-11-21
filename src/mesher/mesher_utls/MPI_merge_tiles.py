import os
import sys
import tempfile
from pathlib import Path
try:
    from mesher.mesher_utls.bootstrap_utils import check_worker_count, ensure_mesher_on_path
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from mesher.mesher_utls.bootstrap_utils import check_worker_count, ensure_mesher_on_path
ensure_mesher_on_path()
import json
import time
import cloudpickle
import numpy as np
from mpi4py import MPI
from mesher.mesher_utls.mpi_progress import complete_task, debug_log


def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True
    return False


def bbox_union(a, b):
    return [
        min(a[0], b[0]),
        min(a[1], b[1]),
        max(a[2], b[2]),
        max(a[3], b[3]),
    ]


def merge_tiles(args):
    if args.get("merge_method", "remesh") != "remesh":
        raise ValueError('Shared-edge meshing has been retired; regenerate tasks with mpi_merge_method="remesh"')
    t0 = time.perf_counter()

    def log_step(msg):
        dt = time.perf_counter() - t0
        debug_log(args, f'[{dt:.2f}s] {msg}')

    tile_a = args['tile_a']
    tile_b = args['tile_b']

    log_step('Loading tile metadata/npz')
    with open(tile_a['meta']) as f:
        meta_a = json.load(f)
    with open(tile_b['meta']) as f:
        meta_b = json.load(f)

    for meta in (meta_a, meta_b):
        if meta.get("merge_method") != "remesh":
            raise ValueError("Cannot merge old shared-edge tiles; regenerate with reuse_mesh=False")

    from mesher.mesher_utls.mpi_remesh import load_domain, load_constraints
    from mesher.mesher_utls.disk_merge import merge_disk_cavity
    bounds_a, bounds_b = meta_a['core_bbox'], meta_b['core_bbox']
    prefix = args['out_prefix']
    archive = Path(prefix + '.npz')
    marker = Path(prefix + '.json')
    with tempfile.TemporaryDirectory(dir=archive.parent, prefix='merge-output-') as directory:
        staged_archive = Path(directory) / archive.name
        staged_marker = Path(directory) / marker.name
        report = merge_disk_cavity(tile_a['npz'], bounds_a, tile_b['npz'], bounds_b,
            load_domain(args['outer_polygon_shp']), load_constraints(args['constraints']),
            args, staged_archive)
        bounds = bbox_union(bounds_a, bounds_b)
        with staged_marker.open('w') as stream:
            json.dump({'core_bbox': bounds, 'tile_bbox': bounds, 'band_width': 0,
                       'empty_tile': report['output_faces'] == 0, 'merge_method': 'remesh',
                       'remesh_report': report}, stream)
        # Invalidate an earlier completion marker only after both new files are ready.
        marker.unlink(missing_ok=True)
        os.replace(staged_archive, archive)
        os.replace(staged_marker, marker)


def main(pickle_file: str, disconnect: bool, expected_workers=None):
    check_worker_count(MPI.COMM_WORLD.size, expected_workers)
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    param_args_split = np.array_split(param_args, MPI.COMM_WORLD.size)

    for args in param_args_split[MPI.COMM_WORLD.rank]:
        merge_tiles(args)
        complete_task(args)

    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
