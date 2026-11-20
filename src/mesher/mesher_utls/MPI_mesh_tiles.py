import os
import sys
try:
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
ensure_mesher_on_path()
import json
import cloudpickle
import numpy as np
from mpi4py import MPI


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
    domain = tile_domain(load_domain(args['outer_polygon_shp']), tile_bbox, band, args['tile_prefix'])
    if domain.is_empty:
        print(f'Tile core owns no domain area; skipping halo meshing [{args["tile_prefix"]}]', flush=True)
        verts, tris = np.empty((0, 3)), np.empty((0, 3), dtype=np.int64)
    else:
        verts, tris = remesh_domain(domain, load_constraints(args['constraints']), args, args['base_dir'])
    prefix = args['base_dir'] + args['tile_prefix']
    np.savez(prefix + '.npz', verts=verts, tris=tris)
    with open(prefix + '.json', 'w') as stream:
        json.dump({'tile_bbox': tile_bbox, 'core_bbox': tile_bbox, 'band_width': band,
                   'empty_tile': len(tris) == 0, 'merge_method': 'remesh'}, stream)


def main(pickle_file: str, disconnect: bool):
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    param_args_split = np.array_split(param_args, MPI.COMM_WORLD.size)

    for args in param_args_split[MPI.COMM_WORLD.rank]:
        mesh_tile(args)

    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
