import sys
import math
import cloudpickle
import numpy as np
from mpi4py import MPI
from osgeo import ogr


def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True
    return False


def compute_tile_grid(nranks):
    best_rows = 1
    best_cols = nranks
    best_diff = abs(best_cols - best_rows)
    for rows in range(1, nranks + 1):
        if nranks % rows != 0:
            continue
        cols = nranks // rows
        diff = abs(cols - rows)
        if diff < best_diff:
            best_rows = rows
            best_cols = cols
            best_diff = diff
    return best_rows, best_cols


def vertex_owner(x, y, xmin, ymin, xmax, ymax, rows, cols):
    width = (xmax - xmin) / cols
    height = (ymax - ymin) / rows
    col = int(math.floor((x - xmin) / width))
    row = int(math.floor((ymax - y) / height))
    if col < 0:
        col = 0
    if row < 0:
        row = 0
    if col >= cols:
        col = cols - 1
    if row >= rows:
        row = rows - 1
    return row * cols + col


def load_polygon_geom(shp_path):
    ds = ogr.Open(shp_path)
    if ds is None:
        raise RuntimeError(f'Unable to open polygon shapefile: {shp_path}')
    layer = ds.GetLayer(0)
    geom_union = None
    for feat in layer:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        geom = geom.Clone()
        geom = geom.MakeValid()
        if geom_union is None:
            geom_union = geom
        else:
            geom_union = geom_union.Union(geom)
    ds = None
    if geom_union is None:
        raise RuntimeError(f'No geometry found in polygon shapefile: {shp_path}')
    geom_union = geom_union.MakeValid()
    geom_union = geom_union.Buffer(0)
    return geom_union


def open_array(path, mmap_mode=None):
    if path.endswith('.npy'):
        return np.load(path, mmap_mode=mmap_mode)
    data = np.load(path)
    if isinstance(data, np.lib.npyio.NpzFile):
        raise RuntimeError('Expected .npy for lloyd input; received .npz')
    return data


def init_output_memmap(out_path, verts):
    out = np.lib.format.open_memmap(out_path, mode='w+', dtype=verts.dtype, shape=verts.shape)
    chunk = 1_000_000
    for start in range(0, len(verts), chunk):
        end = min(start + chunk, len(verts))
        out[start:end] = verts[start:end]
    out.flush()
    return out


def do_lloyd(args):
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    verts_path = args['verts_path']
    tris_path = args['tris_path']
    out_verts_path = args['out_verts_path']
    iterations = int(args.get('iterations', 1))
    outer_polygon_shp = args.get('outer_polygon_shp', None)
    boundary_tol = float(args.get('boundary_tol', 0.0))

    verts = open_array(verts_path, mmap_mode='r')
    tris = open_array(tris_path, mmap_mode='r')

    nverts = len(verts)
    ntris = len(tris)

    v_start = (rank * nverts) // size
    v_end = ((rank + 1) * nverts) // size
    v_slice = verts[v_start:v_end]
    local_min = np.min(v_slice[:, :2], axis=0)
    local_max = np.max(v_slice[:, :2], axis=0)

    xmin = comm.allreduce(local_min[0], op=MPI.MIN)
    ymin = comm.allreduce(local_min[1], op=MPI.MIN)
    xmax = comm.allreduce(local_max[0], op=MPI.MAX)
    ymax = comm.allreduce(local_max[1], op=MPI.MAX)

    rows, cols = compute_tile_grid(size)

    t_start = (rank * ntris) // size
    t_end = ((rank + 1) * ntris) // size

    boundary_geom = None
    if outer_polygon_shp and boundary_tol > 0.0:
        boundary_geom = load_polygon_geom(outer_polygon_shp).GetBoundary()

    if rank == 0:
        init_output_memmap(out_verts_path, verts)
    comm.Barrier()
    out = np.lib.format.open_memmap(out_verts_path, mode='r+', dtype=verts.dtype, shape=verts.shape)

    for _ in range(iterations):
        local_contrib = {}
        for tri in tris[t_start:t_end]:
            v0 = verts[tri[0]]
            v1 = verts[tri[1]]
            v2 = verts[tri[2]]
            cx = (v0[0] + v1[0] + v2[0]) / 3.0
            cy = (v0[1] + v1[1] + v2[1]) / 3.0
            for vid in (tri[0], tri[1], tri[2]):
                if vid not in local_contrib:
                    local_contrib[vid] = [0.0, 0.0, 0]
                local_contrib[vid][0] += cx
                local_contrib[vid][1] += cy
                local_contrib[vid][2] += 1

        send_lists = [[] for _ in range(size)]
        for vid, data in local_contrib.items():
            v = verts[vid]
            owner = vertex_owner(v[0], v[1], xmin, ymin, xmax, ymax, rows, cols)
            send_lists[owner].append((vid, data[0], data[1], data[2]))

        recv_lists = comm.alltoall(send_lists)

        owned_updates = {}
        for entries in recv_lists:
            for vid, sx, sy, count in entries:
                if vid not in owned_updates:
                    owned_updates[vid] = [0.0, 0.0, 0]
                owned_updates[vid][0] += sx
                owned_updates[vid][1] += sy
                owned_updates[vid][2] += count

        for vid, data in owned_updates.items():
            if data[2] <= 0:
                continue
            v = verts[vid]
            if boundary_geom is not None:
                pt = ogr.Geometry(ogr.wkbPoint)
                pt.AddPoint(v[0], v[1])
                if pt.Distance(boundary_geom) <= boundary_tol:
                    continue
            out[vid, 0] = data[0] / data[2]
            out[vid, 1] = data[1] / data[2]
            out[vid, 2] = v[2]
        out.flush()
        comm.Barrier()
        verts = open_array(out_verts_path, mmap_mode='r')
        comm.Barrier()


def main(pickle_file: str, disconnect: bool):
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        args = cloudpickle.load(f)

    do_lloyd(args)

    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
