import os
import sys
import json
import cloudpickle
import subprocess
import numpy as np
from mpi4py import MPI
from osgeo import ogr, gdal, osr


def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True
    return False


def expand_bbox(bbox, buffer_dist):
    return [
        bbox[0] - buffer_dist,
        bbox[1] - buffer_dist,
        bbox[2] + buffer_dist,
        bbox[3] + buffer_dist
    ]


def bbox_intersection(a, b):
    xmin = max(a[0], b[0])
    ymin = max(a[1], b[1])
    xmax = min(a[2], b[2])
    ymax = min(a[3], b[3])
    if xmin >= xmax or ymin >= ymax:
        return None
    return [xmin, ymin, xmax, ymax]


def osr_from_wkt(wkt):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    return srs


def write_polygon_shp(path, geom, srs_wkt):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(path):
        driver.DeleteDataSource(path)

    ds = driver.CreateDataSource(path)
    srs = osr_from_wkt(srs_wkt)
    layer = ds.CreateLayer(path, srs, ogr.wkbPolygon)
    feature_defn = layer.GetLayerDefn()
    feature = ogr.Feature(feature_defn)
    feature.SetGeometry(geom)
    layer.CreateFeature(feature)
    ds = None


def extract_polygon(geom):
    if geom is None:
        return None

    gtype = geom.GetGeometryType()
    if gtype in (ogr.wkbPolygon, ogr.wkbPolygon25D):
        return geom.Clone()
    if gtype in (ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D,
                 ogr.wkbGeometryCollection, ogr.wkbGeometryCollection25D):
        max_area = -1
        poly = None
        for i in range(geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            pg = extract_polygon(g)
            if pg is None:
                continue
            area = pg.GetArea()
            if area > max_area:
                max_area = area
                poly = pg.Clone()
        return poly
    return None


def polygon_exterior_coords(geom):
    poly = extract_polygon(geom)
    if poly is None:
        raise RuntimeError(f'Unsupported geometry type for seam polygon {geom.GetGeometryType()}')

    ring = poly.GetGeometryRef(0)
    coords = []
    for i in range(ring.GetPointCount()):
        x, y, _ = ring.GetPoint(i)
        coords.append([x, y])
    return coords


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


def triangles_to_union_polygon(verts, tris, mask):
    geom_union = None
    for tri in tris[mask]:
        v0 = verts[tri[0]]
        v1 = verts[tri[1]]
        v2 = verts[tri[2]]

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(v0[0], v0[1])
        ring.AddPoint(v1[0], v1[1])
        ring.AddPoint(v2[0], v2[1])
        ring.AddPoint(v0[0], v0[1])

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        if geom_union is None:
            geom_union = poly.Clone()
        else:
            geom_union = geom_union.Union(poly)

    if geom_union is None:
        raise RuntimeError('No seam polygons produced')
    return geom_union


def build_mesher_exec_str(args, poly_file, interior_plgs, points_file):
    execstr = '%s --poly-file %s --tolerance %s --raster %s --area %s --min-area %s --error-metric %s --lloyd %d --interior-plgs-file %s --points-file %s' % \
              (args['mesher_path'],
               poly_file,
               args['max_tolerance'],
               args['dem_path'],
               args['max_area'],
               args['min_area'],
               args['errormetric'],
               args['lloyd_itr'],
               interior_plgs,
               points_file
               )

    if args['is_geographic']:
        execstr += ' --is-geographic true'

    if args['use_weights']:
        execstr += ' --weight %s' % args['topo_weight']
        execstr += ' --weight-threshold %s' % args['weight_threshold']

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


def centroid_mask_in_bbox(verts, tris, bbox):
    tri_pts = verts[tris][:, :, :2]
    centroids = np.mean(tri_pts, axis=1)
    return (
        (centroids[:, 0] >= bbox[0]) &
        (centroids[:, 0] <= bbox[2]) &
        (centroids[:, 1] >= bbox[1]) &
        (centroids[:, 1] <= bbox[3])
    )


def thin_points(points, spacing):
    if spacing <= 0:
        return points

    grid = {}
    kept = []
    inv = 1.0 / spacing
    for p in points:
        key = (int(p[0] * inv), int(p[1] * inv))
        if key in grid:
            continue
        grid[key] = True
        kept.append(p)
    if len(kept) == 0:
        return points
    return np.asarray(kept)


def add_vertex(verts_out, index_map, coord, tol=1e-6):
    key = (round(coord[0] / tol), round(coord[1] / tol))
    if key in index_map:
        return index_map[key]
    idx = len(verts_out)
    verts_out.append([coord[0], coord[1], coord[2]])
    index_map[key] = idx
    return idx


def merge_tiles(args):
    tile_a = args['tile_a']
    tile_b = args['tile_b']

    with open(tile_a['meta']) as f:
        meta_a = json.load(f)
    with open(tile_b['meta']) as f:
        meta_b = json.load(f)

    data_a = np.load(tile_a['npz'])
    data_b = np.load(tile_b['npz'])

    verts_a = data_a['verts']
    tris_a = data_a['tris']
    band_a = data_a['band_tri_mask']
    verts_b = data_b['verts']
    tris_b = data_b['tris']
    band_b = data_b['band_tri_mask']

    band_width = meta_a['band_width']
    bbox_a = meta_a['tile_bbox']
    bbox_b = meta_b['tile_bbox']

    seam_bbox = bbox_intersection(expand_bbox(bbox_a, band_width), expand_bbox(bbox_b, band_width))
    if seam_bbox is None:
        raise RuntimeError('No seam bbox overlap between tiles')

    mask_a = centroid_mask_in_bbox(verts_a, tris_a, seam_bbox)
    mask_b = centroid_mask_in_bbox(verts_b, tris_b, seam_bbox)
    if not np.any(mask_a):
        mask_a = band_a
    if not np.any(mask_b):
        mask_b = band_b

    seam_geom = triangles_to_union_polygon(verts_a, tris_a, mask_a)
    seam_geom = seam_geom.Union(triangles_to_union_polygon(verts_b, tris_b, mask_b))
    seam_poly = extract_polygon(seam_geom)
    if seam_poly is None:
        raise RuntimeError('Unsupported geometry type for seam polygon after union')

    dem_ds = gdal.Open(args['dem_path'])
    if dem_ds is None:
        raise RuntimeError('Unable to open DEM for merge')

    srs_wkt = dem_ds.GetProjection()
    seam_shp = args['out_prefix'] + '_seam.shp'
    write_polygon_shp(seam_shp, seam_poly, srs_wkt)

    coords = polygon_exterior_coords(seam_poly)
    poly_file = args['out_prefix'] + '_seam.poly'
    write_poly_from_coords(poly_file, coords)

    interior_PLGS = {
        "type": "FeatureCollection",
        "name": "interior_PLGS",
        "features": []
    }

    for cpath in args['constraints']:
        outname = args['out_prefix'] + '_constraint_' + os.path.splitext(os.path.basename(cpath))[0]
        exec_str = '%sogr2ogr -f "ESRI Shapefile" -clipsrc %s %s %s' % (
            args['gdal_prefix'], seam_shp, outname + '.shp', cpath)
        subprocess.check_call(exec_str, shell=True)

        exec_str = '%sogr2ogr -f GeoJSON -nlt LINESTRING -explodecollections %s %s' % (
            args['gdal_prefix'], outname + '.geojson', outname + '.shp')
        subprocess.check_call(exec_str, shell=True)

        with open(outname + '.geojson') as f:
            gj = json.load(f)
            for feat in gj.get('features', []):
                if feat.get('geometry') is not None:
                    interior_PLGS['features'].append(feat)

    interior_plgs_file = args['out_prefix'] + '_interior_PLGS.geojson'
    with open(interior_plgs_file, 'w') as fp:
        json.dump(interior_PLGS, fp)

    points_file = args['out_prefix'] + '_points.txt'
    seam_points = []
    for tri in tris_a[mask_a]:
        for idx in tri:
            seam_points.append(verts_a[idx])
    for tri in tris_b[mask_b]:
        for idx in tri:
            seam_points.append(verts_b[idx])

    if len(seam_points) == 0:
        raise RuntimeError('No seam points available for merge')

    seam_points = np.asarray(seam_points)
    seam_spacing = args.get('seam_point_spacing', None)
    if seam_spacing is not None:
        seam_points = thin_points(seam_points, seam_spacing)
    with open(points_file, 'w') as f:
        for v in seam_points:
            f.write(f'{v[0]} {v[1]}\n')

    execstr = build_mesher_exec_str(args, poly_file, interior_plgs_file, points_file)
    subprocess.check_call(execstr, shell=True)

    node_file = poly_file.replace('.poly', '.1.node')
    ele_file = poly_file.replace('.poly', '.1.ele')

    seam_verts = []
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
                    seam_verts.append([float(items[1]), float(items[2]), 0.0])

    seam_tris = []
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
                    seam_tris.append([int(items[1]) - 1, int(items[2]) - 1, int(items[3]) - 1])

    seam_verts = np.asarray(seam_verts, dtype=float)
    seam_tris = np.asarray(seam_tris, dtype=int)

    verts_out = []
    index_map = {}
    tris_out = []

    for tri in tris_a[~mask_a]:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, verts_a[idx]))
        tris_out.append(idxs)

    for tri in tris_b[~mask_b]:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, verts_b[idx]))
        tris_out.append(idxs)

    for tri in seam_tris:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, seam_verts[idx]))
        tris_out.append(idxs)

    verts_out = np.asarray(verts_out, dtype=float)
    tris_out = np.asarray(tris_out, dtype=int)

    merged_bbox = [
        min(bbox_a[0], bbox_b[0]),
        min(bbox_a[1], bbox_b[1]),
        max(bbox_a[2], bbox_b[2]),
        max(bbox_a[3], bbox_b[3])
    ]

    band_mask = compute_band_mask(verts_out, tris_out, merged_bbox, band_width)

    npz_path = args['out_prefix'] + '.npz'
    np.savez(npz_path, verts=verts_out, tris=tris_out, band_tri_mask=band_mask)

    meta_path = args['out_prefix'] + '.json'
    with open(meta_path, 'w') as f:
        json.dump({
            'tile_bbox': merged_bbox,
            'band_width': band_width
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
