import os
import sys
try:
    from mesher.mesher_utls.bootstrap_utils import check_worker_count, ensure_mesher_on_path
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from mesher.mesher_utls.bootstrap_utils import check_worker_count, ensure_mesher_on_path
ensure_mesher_on_path()
from pathlib import Path
import tempfile
import cloudpickle
import numpy as np
from mpi4py import MPI
from mesher.mesher_utls.mpi_resume import publish
from mesher.mesher_utls.mpi_progress import complete_task, debug_log


def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True
    return False


def mesh_tile(args):
    if args.get("merge_method", "remesh") != "remesh":
        raise ValueError('Shared-edge meshing has been retired; regenerate tasks with mpi_merge_method="remesh"')
    from mesher.mesher_utls.mpi_remesh import load_domain, load_constraints, remesh_domain, tile_domain
    tile_bbox = args['tile_bbox']
    band = args['band_width']
    domain = tile_domain(load_domain(args['outer_polygon_shp']), tile_bbox, band, args['tile_prefix'],
                         mesher_debug=args.get('mesher_debug', False))
    if domain.is_empty:
        debug_log(args, f'Tile core owns no domain area; skipping halo meshing [{args["tile_prefix"]}]')
        verts, tris = np.empty((0, 3)), np.empty((0, 3), dtype=np.int64)
    else:
        verts, tris = remesh_domain(domain, load_constraints(args['constraints']), args, args['base_dir'])
    prefix = args['base_dir'] + args['tile_prefix']
    with tempfile.TemporaryDirectory(dir=args['base_dir'], prefix='tile-output-') as work:
        staged = Path(work)/'tile.npz'
        np.savez(staged, verts=verts, tris=tris)
        publish(staged, prefix, dict(tile_bbox=tile_bbox, core_bbox=tile_bbox, band_width=band,
                                     empty_tile=len(tris) == 0, merge_method='remesh'), args)


def main(pickle_file: str, disconnect: bool, expected_workers=None):
    check_worker_count(MPI.COMM_WORLD.size, expected_workers)
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    param_args_split = np.array_split(param_args, MPI.COMM_WORLD.size)

    for args in param_args_split[MPI.COMM_WORLD.rank]:
        mesh_tile(args)
        complete_task(args)

    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
