import os
import sys
import json
import cloudpickle
import subprocess
import numpy as np
from mpi4py import MPI
from osgeo import ogr, gdal, osr

gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()
osr.UseExceptions()


def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True
    return False


def write_bbox_shp(path, bbox, srs_wkt):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(path):
        driver.DeleteDataSource(path)

    ds = driver.CreateDataSource(path)
    srs = osr_from_wkt(srs_wkt)
    layer = ds.CreateLayer(path, srs, ogr.wkbPolygon)
    feature_defn = layer.GetLayerDefn()

    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(bbox[0], bbox[1])
    ring.AddPoint(bbox[2], bbox[1])
    ring.AddPoint(bbox[2], bbox[3])
    ring.AddPoint(bbox[0], bbox[3])
    ring.AddPoint(bbox[0], bbox[1])

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    feature = ogr.Feature(feature_defn)
    feature.SetGeometry(poly)
    layer.CreateFeature(feature)

    ds = None


def osr_from_wkt(wkt):
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    return srs


def longest_linestring_coords(plgs):
    idx = -1
    cmax = -1
    coords_out = None

    for i, features in enumerate(plgs['features']):
        geom = features.get('geometry')
        if geom is None:
            continue

        coords = []
        if geom['type'] == 'LineString':
            coords = geom['coordinates']
        elif geom['type'] == 'MultiLineString':
            len_ml = -1
            idx_ml = -1
            for j, lines in enumerate(geom['coordinates']):
                l = len(lines)
                if l > len_ml:
                    len_ml = l
                    idx_ml = j
            coords = geom['coordinates'][idx_ml]
        else:
            continue

        if len(coords) > cmax:
            cmax = len(coords)
            idx = i
            coords_out = coords

    if idx == -1 or coords_out is None:
        raise RuntimeError('Unable to find a valid linestring for the tile boundary')

    return coords_out


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

    tile_bbox_shp = base_dir + tile_prefix + '_bbox.shp'
    write_bbox_shp(tile_bbox_shp, band_bbox, srs_wkt)

    tile_domain_shp = base_dir + tile_prefix + '_domain.shp'
    exec_str = '%sogr2ogr -f "ESRI Shapefile" -clipsrc %s %s %s' % (
        args['gdal_prefix'], tile_bbox_shp, tile_domain_shp, args['outer_polygon_shp'])
    subprocess.check_call(exec_str, shell=True)

    tile_line_shp = base_dir + tile_prefix + '_line.shp'
    exec_str = '%sogr2ogr -overwrite %s %s -nlt LINESTRING' % (
        args['gdal_prefix'], tile_line_shp, tile_domain_shp)
    subprocess.check_call(exec_str, shell=True)

    tile_geojson = base_dir + tile_prefix + '_boundary.geojson'
    exec_str = '%sogr2ogr -f GeoJSON %s %s' % (
        args['gdal_prefix'], tile_geojson, tile_line_shp)
    subprocess.check_call(exec_str, shell=True)

    with open(tile_geojson) as f:
        plgs = json.load(f)

    coords = longest_linestring_coords(plgs)

    poly_file = base_dir + tile_prefix + '.poly'
    write_poly_from_coords(poly_file, coords)

    interior_PLGS = {
        "type": "FeatureCollection",
        "name": "interior_PLGS",
        "features": []
    }

    for cpath in args['constraints']:
        outname = base_dir + tile_prefix + '_constraint_' + os.path.splitext(os.path.basename(cpath))[0]
        exec_str = '%sogr2ogr -f "ESRI Shapefile" -clipsrc %s %s %s' % (
            args['gdal_prefix'], tile_bbox_shp, outname + '.shp', cpath)
        subprocess.check_call(exec_str, shell=True)

        exec_str = '%sogr2ogr -f GeoJSON -nlt LINESTRING -explodecollections %s %s' % (
            args['gdal_prefix'], outname + '.geojson', outname + '.shp')
        subprocess.check_call(exec_str, shell=True)

        with open(outname + '.geojson') as f:
            gj = json.load(f)
            for feat in gj.get('features', []):
                if feat.get('geometry') is not None:
                    interior_PLGS['features'].append(feat)

    interior_plgs_file = base_dir + tile_prefix + '_interior_PLGS.geojson'
    with open(interior_plgs_file, 'w') as fp:
        json.dump(interior_PLGS, fp)

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
gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()
osr.UseExceptions()
