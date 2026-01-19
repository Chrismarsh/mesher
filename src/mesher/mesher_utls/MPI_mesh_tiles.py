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
import subprocess
import numpy as np
from mpi4py import MPI
from osgeo import ogr, gdal
from mesher.mesher_utls.ogr_utils import bbox_to_polygon_geom, load_polygon_geom, normalize_polygon, \
    linestring_features_from_geom

gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()


def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True
    return False


def longest_linestring_coords_from_geom(geom):
    lines = linestring_features_from_geom(geom)
    coords_out = None
    cmax = -1
    for feat in lines:
        coords = feat["geometry"]["coordinates"]
        if len(coords) > cmax:
            cmax = len(coords)
            coords_out = coords
    if coords_out is None:
        coords_out = polygon_exterior_coords_from_geom(geom)
    if coords_out is None:
        gtype = geom.GetGeometryType() if geom is not None else None
        gname = ogr.GeometryTypeToName(gtype) if gtype is not None else 'None'
        raise RuntimeError(
            f'Unable to find a valid linestring for the tile boundary '
            f'(geom_type={gname})'
        )
    return coords_out


def polygon_exterior_coords_from_geom(geom):
    if geom is None:
        return None
    if geom.IsEmpty():
        return None
    gtype = geom.GetGeometryType()
    if gtype in (ogr.wkbPolygon, ogr.wkbPolygon25D):
        ring = geom.GetGeometryRef(0)
        if ring is not None and ring.GetPointCount() > 0:
            return ring.GetPoints()
        env = geom.GetEnvelope()
        if env[0] < env[1] and env[2] < env[3]:
            return [
                (env[0], env[2]),
                (env[1], env[2]),
                (env[1], env[3]),
                (env[0], env[3]),
                (env[0], env[2]),
            ]
    if gtype in (ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D):
        best = None
        best_area = -1.0
        for i in range(geom.GetGeometryCount()):
            poly = geom.GetGeometryRef(i)
            if poly is None:
                continue
            area = poly.GetArea()
            if area > best_area:
                best_area = area
                best = poly
        if best is not None:
            ring = best.GetGeometryRef(0)
            if ring is not None and ring.GetPointCount() > 0:
                return ring.GetPoints()
            env = best.GetEnvelope()
            if env[0] < env[1] and env[2] < env[3]:
                return [
                    (env[0], env[2]),
                    (env[1], env[2]),
                    (env[1], env[3]),
                    (env[0], env[3]),
                    (env[0], env[2]),
                ]
    if gtype == ogr.wkbLinearRing:
        if geom.GetPointCount() > 0:
            return geom.GetPoints()
    return None


def write_poly_from_coords(poly_path, coords):
    with open(poly_path, 'w') as f:
        header = '%d 2 0 0\n' % (len(coords))
        f.write(header)
        vert = 1
        for c in coords:
            f.write('%d %17.11f %17.11f\n' % (vert, c[0], c[1]))
            vert += 1

        f.write('\n')
        header = '%d 0\n' % (len(coords))
        f.write(header)

        for i in range(len(coords)):
            if i + 1 == len(coords):
                f.write('%d %d %d\n' % (i + 1, i + 1, 1))
            else:
                f.write('%d %d %d\n' % (i + 1, i + 1, i + 2))

        f.write('0\n')


def build_mesher_exec_str(args, poly_file, interior_plgs):
    execstr = '%s --poly-file %s --tolerance %s --raster %s --area %s --min-area %s --error-metric %s --lloyd %d --interior-plgs-file %s' % \
              (args['mesher_path'],
               poly_file,
               args['max_tolerance'],
               args['dem_path'],
               args['max_area'],
               args['min_area'],
               args['errormetric'],
               args['lloyd_itr'],
               interior_plgs
               )

    if args['is_geographic']:
        execstr += ' --is-geographic true'

    if args['use_weights']:
        execstr += ' --weight %s' % args['topo_weight']
        execstr += ' --weight-threshold %s' % args['weight_threshold']
    if args.get('mesher_debug', False):
        execstr += ' --debug true'

    for key, data in args['parameter_files'].items():
        if 'tolerance' in data:
            if data['method'] == 'mode':
                execstr += ' --category-raster %s --category-frac %s' % (data['filename'][0], data['tolerance'])
            else:
                execstr += ' --raster %s --tolerance %s' % (data['filename'][0], data['tolerance'])
        if args['use_weights'] and 'weight' in data:
            execstr += ' --weight %s' % data['weight']

    for key, data in args['initial_conditions'].items():
        if 'tolerance' in data:
            if data['method'] == 'mode':
                execstr += ' --category-raster %s --category-frac %s' % (data['filename'], data['tolerance'])
            else:
                execstr += ' --raster %s --tolerance %s' % (data['filename'][0], data['tolerance'])
        if args['use_weights'] and 'weight' in data:
            execstr += ' --weight %s' % data['weight']

    return execstr


def compute_band_mask(verts, tris, tile_bbox, band_width):
    interior_bbox = [
        tile_bbox[0] + band_width,
        tile_bbox[1] + band_width,
        tile_bbox[2] - band_width,
        tile_bbox[3] - band_width
    ]

    if interior_bbox[0] >= interior_bbox[2] or interior_bbox[1] >= interior_bbox[3]:
        return np.ones(len(tris), dtype=bool)

    tri_pts = verts[tris][:, :, :2]
    centroids = np.mean(tri_pts, axis=1)

    inside = (
        (centroids[:, 0] >= interior_bbox[0]) &
        (centroids[:, 0] <= interior_bbox[2]) &
        (centroids[:, 1] >= interior_bbox[1]) &
        (centroids[:, 1] <= interior_bbox[3])
    )
    return ~inside


def mesh_tile(args):
    tile_prefix = args['tile_prefix']
    base_dir = args['base_dir']
    tile_bbox = args['tile_bbox']
    band_width = args['band_width']

    band_bbox = [
        tile_bbox[0] - band_width,
        tile_bbox[1] - band_width,
        tile_bbox[2] + band_width,
        tile_bbox[3] + band_width
    ]

    dem_ds = gdal.Open(args['dem_path'])
    if dem_ds is None:
        raise RuntimeError('Unable to open DEM for tile meshing')

    srs_wkt = dem_ds.GetProjection()

    bbox_poly = bbox_to_polygon_geom(band_bbox)
    outer_poly = load_polygon_geom(args['outer_polygon_shp'])
    tile_poly = normalize_polygon(outer_poly.Intersection(bbox_poly))
    if tile_poly is None:
        raise RuntimeError('Tile polygon invalid after bbox clip')
    if tile_poly.IsEmpty():
        verts = np.zeros((0, 3), dtype=float)
        tris = np.zeros((0, 3), dtype=int)
        band_mask = np.zeros((0,), dtype=bool)
        npz_path = base_dir + tile_prefix + '.npz'
        np.savez(npz_path, verts=verts, tris=tris, band_tri_mask=band_mask)
        meta_path = base_dir + tile_prefix + '.json'
        with open(meta_path, 'w') as f:
            json.dump({
                'tile_bbox': tile_bbox,
                'core_bbox': tile_bbox,
                'band_width': band_width,
                'empty_tile': True
            }, f)
        print(f'Tile polygon empty after bbox clip; wrote empty tile {tile_prefix}')
        return
    coords = longest_linestring_coords_from_geom(tile_poly)

    poly_file = base_dir + tile_prefix + '.poly'
    write_poly_from_coords(poly_file, coords)

    interior_PLGS = {
        "type": "FeatureCollection",
        "name": "interior_PLGS",
        "features": []
    }

    for cpath in args['constraints']:
        ds = ogr.Open(cpath)
        if ds is None:
            raise RuntimeError(f'Unable to open constraint {cpath}')
        layer = ds.GetLayer(0)
        for feat in layer:
            geom = feat.GetGeometryRef()
            if geom is None:
                continue
            clipped = geom.Intersection(bbox_poly)
            if clipped is None:
                continue
            interior_PLGS['features'].extend(linestring_features_from_geom(clipped))
        ds = None

    interior_plgs_file = base_dir + tile_prefix + '_interior_PLGS.geojson'
    with open(interior_plgs_file, 'w') as fp:
        json.dump(interior_PLGS, fp)

    if args.get('dump_poly_only', False):
        print(f'Dumping tile poly only (no mesher run): {poly_file}')
        return
    if args.get('dump_poly_files', False):
        print(f'Dumping tile poly files (mesher will still run): {poly_file}')

    execstr = build_mesher_exec_str(args, poly_file, interior_plgs_file)
    subprocess.check_call(execstr, shell=True)

    node_file = poly_file.replace('.poly', '.1.node')
    ele_file = poly_file.replace('.poly', '.1.ele')

    verts = []
    with open(node_file) as f:
        read_header = False
        for line in f:
            if '#' not in line:
                if not read_header:
                    read_header = True
                else:
                    items = line.split()
                    if len(items) < 3:
                        continue
                    verts.append([float(items[1]), float(items[2]), 0.0])

    tris = []
    with open(ele_file) as f:
        read_header = False
        for line in f:
            if '#' not in line:
                if not read_header:
                    read_header = True
                else:
                    items = line.split()
                    if len(items) < 4:
                        continue
                    tris.append([int(items[1]) - 1, int(items[2]) - 1, int(items[3]) - 1])

    verts = np.asarray(verts, dtype=float)
    tris = np.asarray(tris, dtype=int)
    band_mask = compute_band_mask(verts, tris, tile_bbox, band_width)

    npz_path = base_dir + tile_prefix + '.npz'
    np.savez(npz_path, verts=verts, tris=tris, band_tri_mask=band_mask)

    meta_path = base_dir + tile_prefix + '.json'
    with open(meta_path, 'w') as f:
        json.dump({
            'tile_bbox': tile_bbox,
            'core_bbox': tile_bbox,
            'band_width': band_width
        }, f)


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
