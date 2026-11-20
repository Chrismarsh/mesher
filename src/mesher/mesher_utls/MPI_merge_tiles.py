import os
import sys
try:
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
ensure_mesher_on_path()
import json
import math
import shutil
import time
import cloudpickle
import numpy as np
from mpi4py import MPI


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


def add_vertex(verts_out, index_map, spatial_map, coord, tol=1e-6):
    key = (round(coord[0] / tol), round(coord[1] / tol))
    if key in index_map:
        return index_map[key]

    if spatial_map is not None and len(verts_out) > 0:
        cell = tol
        cx = int(math.floor(coord[0] / cell))
        cy = int(math.floor(coord[1] / cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for idx in spatial_map.get((cx + dx, cy + dy), []):
                    v = verts_out[idx]
                    if (v[0] - coord[0]) ** 2 + (v[1] - coord[1]) ** 2 <= tol * tol:
                        return idx

    idx = len(verts_out)
    verts_out.append([coord[0], coord[1], coord[2]])
    index_map[key] = idx
    if spatial_map is not None:
        cell = tol
        cx = int(math.floor(coord[0] / cell))
        cy = int(math.floor(coord[1] / cell))
        spatial_map.setdefault((cx, cy), []).append(idx)
    return idx


def merge_tiles(args):
    t0 = time.perf_counter()

    def log_step(msg):
        dt = time.perf_counter() - t0
        print(f'[{dt:.2f}s] {msg}', flush=True)

    tile_a = args['tile_a']
    tile_b = args['tile_b']

    log_step('Loading tile metadata/npz')
    with open(tile_a['meta']) as f:
        meta_a = json.load(f)
    with open(tile_b['meta']) as f:
        meta_b = json.load(f)

    if args.get('merge_method') == 'remesh':
        from mesher.mesher_utls.mpi_remesh import load_domain, load_constraints, merge_cavity
        import tempfile
        bounds_a, bounds_b = meta_a['core_bbox'], meta_b['core_bbox']
        with np.load(tile_a['npz']) as a, np.load(tile_b['npz']) as b:
            av, at, bv, bt = a['verts'], a['tris'], b['verts'], b['tris']
            input_seconds = time.perf_counter() - t0
            verts, tris, report = merge_cavity(av, at, bounds_a,
                bv, bt, bounds_b, load_domain(args['outer_polygon_shp']),
                load_constraints(args['constraints']), args, os.path.dirname(args['out_prefix']) or '.')
        # All checks completed before publishing either file. Metadata is the
        # completion marker and is written last; failed candidates stay temporary.
        report['input_seconds'] = input_seconds
        write_started = time.perf_counter()
        prefix = args['out_prefix']
        with tempfile.NamedTemporaryFile(dir=os.path.dirname(prefix) or '.', suffix='.npz', delete=False) as stream:
            temporary = stream.name
        try:
            np.savez(temporary, verts=verts, tris=tris)
            os.replace(temporary, prefix + '.npz')
        finally:
            if os.path.exists(temporary): os.unlink(temporary)
        report['output_npz_seconds'] = time.perf_counter() - write_started
        with open(prefix + '.json.tmp', 'w') as stream:
            json.dump({'core_bbox': bbox_union(bounds_a, bounds_b),
                       'tile_bbox': bbox_union(bounds_a, bounds_b), 'band_width': 0,
                       'empty_tile': len(tris) == 0, 'merge_method': 'remesh', 'remesh_report': report}, stream)
        os.replace(prefix + '.json.tmp', prefix + '.json')
        return

    if meta_a.get('empty_tile') or meta_b.get('empty_tile'):
        out_prefix = args['out_prefix']
        out_npz = out_prefix + '.npz'
        out_meta = out_prefix + '.json'
        if meta_a.get('empty_tile') and meta_b.get('empty_tile'):
            verts = np.zeros((0, 3), dtype=float)
            tris = np.zeros((0, 3), dtype=int)
            np.savez(out_npz, verts=verts, tris=tris)
            with open(out_meta, 'w') as f:
                json.dump({
                    'tile_bbox': meta_a.get('tile_bbox', meta_b.get('tile_bbox')),
                    'core_bbox': meta_a.get('core_bbox', meta_b.get('core_bbox')),
                    'band_width': 0.0,
                    'empty_tile': True
                }, f)
            print(f'Merge skipped: both tiles empty for {out_prefix}')
            return
        if meta_a.get('empty_tile'):
            shutil.copyfile(tile_b['npz'], out_npz)
            shutil.copyfile(tile_b['meta'], out_meta)
            print(f'Merge skipped: tile_a empty for {out_prefix}')
            return
        shutil.copyfile(tile_a['npz'], out_npz)
        shutil.copyfile(tile_a['meta'], out_meta)
        print(f'Merge skipped: tile_b empty for {out_prefix}')
        return

    data_a = np.load(tile_a['npz'])
    data_b = np.load(tile_b['npz'])

    log_step('Merging tiles')
    verts_a = data_a['verts']
    tris_a = data_a['tris']
    verts_b = data_b['verts']
    tris_b = data_b['tris']
    snap_tol = args.get('merge_snap_tol', 1e-6) or 1e-6

    verts_out = []
    index_map = {}
    spatial_map = {}
    tris_out = []

    for tri in tris_a:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, spatial_map, verts_a[idx], snap_tol))
        tris_out.append(idxs)

    for tri in tris_b:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, spatial_map, verts_b[idx], snap_tol))
        tris_out.append(idxs)

    verts_out = np.asarray(verts_out, dtype=float)
    tris_out = np.asarray(tris_out, dtype=int)

    out_prefix = args['out_prefix']
    npz_path = out_prefix + '.npz'
    np.savez(npz_path, verts=verts_out, tris=tris_out)

    meta_path = out_prefix + '.json'
    bbox_a = meta_a.get('tile_bbox', meta_a.get('core_bbox'))
    bbox_b = meta_b.get('tile_bbox', meta_b.get('core_bbox'))
    with open(meta_path, 'w') as f:
        json.dump({
            'tile_bbox': bbox_union(bbox_a, bbox_b),
            'core_bbox': bbox_union(bbox_a, bbox_b),
            'band_width': 0.0
        }, f)


def main(pickle_file: str, disconnect: bool):
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    param_args_split = np.array_split(param_args, MPI.COMM_WORLD.size)

    for args in param_args_split[MPI.COMM_WORLD.rank]:
        merge_tiles(args)

    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
