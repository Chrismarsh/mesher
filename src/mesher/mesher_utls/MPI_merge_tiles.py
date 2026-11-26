import os
import hashlib
from contextlib import suppress
import shutil
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
from mesher.mesher_utls.mpi_resume import publish
from mesher.mesher_utls.mpi_progress import complete_task, debug_log, report_phase


# A site launcher passes the disconnect flag as text. Only an explicit true means this worker
# has a spawned MPI parent.
def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True
    return False


# The next merge owns the combined rectangle of both cores, including portions with no
# triangles.
def bbox_union(a, b):
    return [
        min(a[0], b[0]),
        min(a[1], b[1]),
        max(a[2], b[2]),
        max(a[3], b[3]),
    ]


# Read the parents' committed core bounds, solve their common seam, and publish the staged mesh
# and marker together for the next merge level.
def merge_tiles(args):
    if args.get("merge_method", "remesh") != "remesh":
        raise ValueError('Shared-edge meshing has been retired; regenerate tasks with mpi_merge_method="remesh"')
    t0 = time.perf_counter()

    # Measure each diagnostic from the start of this pair so slow stages can be located even
    # when workers' logs interleave.
    def log_step(msg):
        dt = time.perf_counter() - t0
        debug_log(args, f'[{dt:.2f}s] {msg}')

    tile_a = args['tile_a']
    tile_b = args['tile_b']

    log_step('Loading tile metadata/npz')
    report_phase(args, 'loading merge metadata and constraints')
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
    with tempfile.TemporaryDirectory(dir=archive.parent, prefix='merge-output-') as directory:
        staged_archive = Path(directory) / archive.name
        cache = None
        if '_cavity_team' in args:
            # A fresh run, a regenerated parent or a rebuilt native solver must
            # never pick up pieces produced under the old criteria or geometry.
            cache = Path(directory) / 'pieces'
            if 'checkpoint_run' in args:
                binary = Path(shutil.which(args['mesher_path']) or args['mesher_path']).resolve()
                stat = binary.stat()
                identity = [1, args['checkpoint_run'], meta_a['checkpoint']['generation'],
                            meta_b['checkpoint']['generation'], str(binary), stat.st_size, stat.st_mtime_ns,
                            args.get('remesh_max_vertices', 500000), args.get('remesh_max_active_faces', 2000000),
                            args.get('remesh_max_attempts', 12), args['remesh_band_width']]
                key = hashlib.sha256(json.dumps(identity).encode()).hexdigest()
                cache = Path(prefix + '.cavity-pieces') / key
                args = dict(args, _seam_retry_path=cache / 'seam-retry.npz')
            args = dict(args, _cavity_cache=cache)

        report = merge_disk_cavity(tile_a['npz'], bounds_a, tile_b['npz'], bounds_b,
            load_domain(args['outer_polygon_shp']), load_constraints(args['constraints']),
            args, staged_archive)
        bounds = bbox_union(bounds_a, bounds_b)
        report_phase(args, 'publishing the merge checkpoint')
        publish(staged_archive, prefix, dict(core_bbox=bounds, tile_bbox=bounds, band_width=0,
                empty_tile=report['output_faces'] == 0, merge_method='remesh', remesh_report=report), args)
        # Until the pair marker is committed a PBS timeout must leave all
        # accepted pieces available. Reproduction logs live outside the cache.
        if cache is not None and cache.exists():
            shutil.rmtree(cache, ignore_errors=True)
            if 'checkpoint_run' in args:
                with suppress(OSError):
                    cache.parent.rmdir()  # Keep any other generation still present.


# Each rank takes a slice of the pending pairs. Mark a pair complete only after its archive and
# metadata have been published.
def main(pickle_file: str, disconnect: bool, expected_workers=None):
    check_worker_count(MPI.COMM_WORLD.size, expected_workers)
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    world = MPI.COMM_WORLD
    groups = min(len(param_args), world.size)
    cap = param_args[0].get('cavity_workers', 1)
    if cap == 1:
        # This also covers fixed-size launchers with more ranks than tasks.
        for args in np.array_split(param_args, world.size)[world.rank]:
            merge_tiles(args)
            complete_task(args)
    else:
        from mesher.mesher_utls.parallel_cavity import CavityTeam
        group = world.rank % groups
        active = world.rank < groups * cap
        comm = world.Split(group if active else MPI.UNDEFINED, world.rank)
        if active:
            try:
                for args in param_args[group::groups]:
                    team = CavityTeam(comm, args)
                    if comm.rank == 0:
                        try:
                            merge_tiles(dict(args, _cavity_team=team))
                            complete_task(args)
                        finally:
                            team.stop()
                    else:
                        team.serve()
            finally:
                comm.Free()

    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
