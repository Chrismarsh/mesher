import os
import sys
import json
import math
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


def extract_polygons(geom):
    if geom is None:
        return []

    gtype = geom.GetGeometryType()
    if gtype in (ogr.wkbPolygon, ogr.wkbPolygon25D):
        return [geom.Clone()]
    if gtype in (ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D,
                 ogr.wkbGeometryCollection, ogr.wkbGeometryCollection25D):
        polys = []
        for i in range(geom.GetGeometryCount()):
            g = geom.GetGeometryRef(i)
            polys.extend(extract_polygons(g))
        return polys
    return []


def clean_ring_coords(coords, tol):
    if len(coords) < 4:
        return coords

    cleaned = [coords[0]]
    for pt in coords[1:]:
        if abs(pt[0] - cleaned[-1][0]) > tol or abs(pt[1] - cleaned[-1][1]) > tol:
            cleaned.append(pt)

    if len(cleaned) < 4:
        return cleaned

    def collinear(a, b, c):
        area = abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        return area <= tol * tol

    simplified = [cleaned[0]]
    for i in range(1, len(cleaned) - 1):
        if not collinear(cleaned[i - 1], cleaned[i], cleaned[i + 1]):
            simplified.append(cleaned[i])
    simplified.append(cleaned[-1])

    if simplified[0] != simplified[-1]:
        simplified.append(simplified[0])
    return simplified


def polygon_exterior_coords(geom):
    polys = extract_polygons(geom)
    if not polys:
        raise RuntimeError(f'Unsupported geometry type for seam polygon {geom.GetGeometryType()}')
    poly = polys[0]

    ring = poly.GetGeometryRef(0)
    coords = []
    for i in range(ring.GetPointCount()):
        x, y, _ = ring.GetPoint(i)
        coords.append([x, y])
    return coords


def bbox_to_polygon(bbox):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(bbox[0], bbox[1])
    ring.AddPoint(bbox[2], bbox[1])
    ring.AddPoint(bbox[2], bbox[3])
    ring.AddPoint(bbox[0], bbox[3])
    ring.AddPoint(bbox[0], bbox[1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


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


def centroid_mask_outside_bbox(verts, tris, bbox):
    return ~centroid_mask_in_bbox(verts, tris, bbox)


def tri_inside_bbox(verts, tris, bbox):
    v = verts[tris][:, :, :2]
    inside = (
        (v[:, :, 0] >= bbox[0]) &
        (v[:, :, 0] <= bbox[2]) &
        (v[:, :, 1] >= bbox[1]) &
        (v[:, :, 1] <= bbox[3])
    )
    return np.all(inside, axis=1)


def dedupe_points(points, tol):
    if len(points) == 0:
        return points
    grid = {}
    kept = []
    inv = 1.0 / tol if tol > 0 else 1.0
    for p in points:
        key = (round(p[0] * inv), round(p[1] * inv))
        if key in grid:
            continue
        grid[key] = True
        kept.append(p)
    return np.asarray(kept)


def shared_core_edge(bbox_a, bbox_b, eps=1e-6):
    # Returns (orientation, const, minv, maxv) or None.
    # orientation: "v" for x=const with y in [minv, maxv], "h" for y=const with x in [minv, maxv]
    if abs(bbox_a[2] - bbox_b[0]) <= eps:
        y0 = max(bbox_a[1], bbox_b[1])
        y1 = min(bbox_a[3], bbox_b[3])
        if y0 < y1:
            return ("v", bbox_a[2], y0, y1)
    if abs(bbox_b[2] - bbox_a[0]) <= eps:
        y0 = max(bbox_a[1], bbox_b[1])
        y1 = min(bbox_a[3], bbox_b[3])
        if y0 < y1:
            return ("v", bbox_a[0], y0, y1)
    if abs(bbox_a[3] - bbox_b[1]) <= eps:
        x0 = max(bbox_a[0], bbox_b[0])
        x1 = min(bbox_a[2], bbox_b[2])
        if x0 < x1:
            return ("h", bbox_a[3], x0, x1)
    if abs(bbox_b[3] - bbox_a[1]) <= eps:
        x0 = max(bbox_a[0], bbox_b[0])
        x1 = min(bbox_a[2], bbox_b[2])
        if x0 < x1:
            return ("h", bbox_a[1], x0, x1)
    return None


def triangle_intersects_polygon(verts, tris, poly):
    mask = np.zeros(len(tris), dtype=bool)
    if len(tris) == 0:
        return mask

    poly_env = poly.GetEnvelope()
    env = (poly_env[0], poly_env[1], poly_env[2], poly_env[3])

    for i, tri in enumerate(tris):
        v0 = verts[tri[0]]
        v1 = verts[tri[1]]
        v2 = verts[tri[2]]

        tri_xmin = min(v0[0], v1[0], v2[0])
        tri_xmax = max(v0[0], v1[0], v2[0])
        tri_ymin = min(v0[1], v1[1], v2[1])
        tri_ymax = max(v0[1], v1[1], v2[1])

        if tri_xmax < env[0] or tri_xmin > env[1] or tri_ymax < env[2] or tri_ymin > env[3]:
            continue

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(v0[0], v0[1])
        ring.AddPoint(v1[0], v1[1])
        ring.AddPoint(v2[0], v2[1])
        ring.AddPoint(v0[0], v0[1])
        poly_tri = ogr.Geometry(ogr.wkbPolygon)
        poly_tri.AddGeometry(ring)

        if poly_tri.Intersects(poly):
            mask[i] = True

    return mask


def thin_points(points, spacing):
    if spacing <= 0:
        return points

    grid = {}
    kept = []
    inv = 1.0 / spacing
    for p in points:
        key = (math.floor(p[0] * inv), math.floor(p[1] * inv))
        if key in grid:
            continue
        grid[key] = True
        kept.append(p)
    if len(kept) == 0:
        return points
    return np.asarray(kept)


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
    bbox_a = meta_a.get('core_bbox', meta_a['tile_bbox'])
    bbox_b = meta_b.get('core_bbox', meta_b['tile_bbox'])

    seam_bbox = bbox_intersection(expand_bbox(bbox_a, band_width), expand_bbox(bbox_b, band_width))
    if seam_bbox is None:
        raise RuntimeError('No seam bbox overlap between tiles')

    shared_edge = shared_core_edge(bbox_a, bbox_b)
    if shared_edge is not None:
        orient, const, v0, v1 = shared_edge
        if orient == "v":
            seam_strip = bbox_to_polygon([const - band_width, v0, const + band_width, v1])
        else:
            seam_strip = bbox_to_polygon([v0, const - band_width, v1, const + band_width])
    else:
        seam_strip = bbox_to_polygon(seam_bbox)

    dem_ds = gdal.Open(args['dem_path'])
    if dem_ds is None:
        raise RuntimeError('Unable to open DEM for merge')

    gt = dem_ds.GetGeoTransform()
    raster_xmax = gt[0] + gt[1] * dem_ds.RasterXSize
    raster_ymin = gt[3] + gt[5] * dem_ds.RasterYSize
    raster_bounds = [gt[0], raster_ymin, raster_xmax, gt[3]]

    raster_ring = ogr.Geometry(ogr.wkbLinearRing)
    raster_ring.AddPoint(raster_bounds[0], raster_bounds[1])
    raster_ring.AddPoint(raster_bounds[2], raster_bounds[1])
    raster_ring.AddPoint(raster_bounds[2], raster_bounds[3])
    raster_ring.AddPoint(raster_bounds[0], raster_bounds[3])
    raster_ring.AddPoint(raster_bounds[0], raster_bounds[1])
    raster_poly = ogr.Geometry(ogr.wkbPolygon)
    raster_poly.AddGeometry(raster_ring)

    seam_strip = seam_strip.Intersection(raster_poly)
    seam_strip = seam_strip.MakeValid()
    seam_strip = seam_strip.Buffer(0)
    seam_strip = extract_polygons(seam_strip)
    seam_strip = seam_strip[0] if seam_strip else None
    if seam_strip is None:
        raise RuntimeError('Seam strip clipped outside raster bounds')

    # remove any triangles from either tile that are in the seam strip or outside ownership
    # ownership is centroid-based to avoid dropping triangles that straddle core bounds
    inside_a = centroid_mask_in_bbox(verts_a, tris_a, bbox_a)
    inside_b = centroid_mask_in_bbox(verts_b, tris_b, bbox_b)

    # initial seam selection based on strip around shared edge
    seam_select_a = triangle_intersects_polygon(verts_a, tris_a, seam_strip)
    seam_select_b = triangle_intersects_polygon(verts_b, tris_b, seam_strip)

    # removal mask includes seam strip and non-owned triangles
    mask_a = seam_select_a | (~inside_a)
    mask_b = seam_select_b | (~inside_b)
    if not np.any(mask_a):
        mask_a = band_a
    if not np.any(mask_b):
        mask_b = band_b

    srs_wkt = dem_ds.GetProjection()

    interior_PLGS = {
        "type": "FeatureCollection",
        "name": "interior_PLGS",
        "features": []
    }

    if shared_edge is not None:
        orient, const, v0, v1 = shared_edge
        if orient == "v":
            coords = [[const, v0], [const, v1]]
        else:
            coords = [[v0, const], [v1, const]]
        interior_PLGS['features'].append({
            "type": "Feature",
            "properties": {"name": "core_shared_edge"},
            "geometry": {"type": "LineString", "coordinates": coords}
        })

    if args.get('seam_constraints', True):
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

    snap_tol = args.get('merge_snap_tol', None)
    seam_spacing = args.get('seam_point_spacing', None)
    if snap_tol is None:
        if seam_spacing is None:
            snap_tol = 1e-6
        else:
            snap_tol = max(1e-6, seam_spacing * 0.01)

    # Build an irregular seam polygon from the triangles we will remove.
    if np.any(mask_a) or np.any(mask_b):
        seam_geom = triangles_to_union_polygon(verts_a, tris_a, mask_a)
        seam_geom = seam_geom.Union(triangles_to_union_polygon(verts_b, tris_b, mask_b))
        seam_poly = extract_polygons(seam_geom)
        seam_poly = seam_poly[0] if seam_poly else None
    else:
        seam_poly = seam_strip

    seam_poly = seam_poly.Intersection(raster_poly)
    seam_poly = seam_poly.MakeValid()
    seam_poly = seam_poly.Buffer(0)
    seam_poly = extract_polygons(seam_poly)
    seam_poly = seam_poly[0] if seam_poly else None
    if seam_poly is None:
        raise RuntimeError('Seam polygon invalid after bounds clipping')

    # Seam selection matches removal mask to keep fill and removal aligned.
    seam_select_a = mask_a
    seam_select_b = mask_b

    seam_shp = args['out_prefix'] + '_seam.shp'
    write_polygon_shp(seam_shp, seam_poly, srs_wkt)

    if args.get('debug_seam', False):
        driver = ogr.GetDriverByName('ESRI Shapefile')
        debug_rm_a = args['out_prefix'] + '_removed_a.shp'
        if os.path.exists(debug_rm_a):
            driver.DeleteDataSource(debug_rm_a)
        ds_a = driver.CreateDataSource(debug_rm_a)
        layer_a = ds_a.CreateLayer('removed_a', osr_from_wkt(srs_wkt), ogr.wkbPolygon)

        debug_rm_b = args['out_prefix'] + '_removed_b.shp'
        if os.path.exists(debug_rm_b):
            driver.DeleteDataSource(debug_rm_b)
        ds_b = driver.CreateDataSource(debug_rm_b)
        layer_b = ds_b.CreateLayer('removed_b', osr_from_wkt(srs_wkt), ogr.wkbPolygon)

        def add_tri(layer, v0, v1, v2):
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(v0[0], v0[1])
            ring.AddPoint(v1[0], v1[1])
            ring.AddPoint(v2[0], v2[1])
            ring.AddPoint(v0[0], v0[1])
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            feat = ogr.Feature(layer.GetLayerDefn())
            feat.SetGeometry(poly)
            layer.CreateFeature(feat)

        for tri in tris_a[mask_a]:
            v0, v1, v2 = verts_a[tri[0]], verts_a[tri[1]], verts_a[tri[2]]
            add_tri(layer_a, v0, v1, v2)

        for tri in tris_b[mask_b]:
            v0, v1, v2 = verts_b[tri[0]], verts_b[tri[1]], verts_b[tri[2]]
            add_tri(layer_b, v0, v1, v2)

        ds_a = None
        ds_b = None

    coords = polygon_exterior_coords(seam_poly)
    seam_simplify_tol = args.get('seam_simplify_tol', None)
    if seam_simplify_tol is not None:
        coords = clean_ring_coords(coords, float(seam_simplify_tol))
        if len(coords) < 4:
            raise RuntimeError('Seam polygon simplified to too few points')
    poly_file = args['out_prefix'] + '_seam.poly'
    write_poly_from_coords(poly_file, coords)

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
    if seam_spacing is not None:
        before = len(seam_points)
        seam_points = thin_points(seam_points, seam_spacing)
        after = len(seam_points)
        print(f'Seam thinning spacing={seam_spacing} kept {after}/{before} points')

    in_bounds = (
        (seam_points[:, 0] >= raster_bounds[0]) &
        (seam_points[:, 0] <= raster_bounds[2]) &
        (seam_points[:, 1] >= raster_bounds[1]) &
        (seam_points[:, 1] <= raster_bounds[3])
    )
    seam_points = seam_points[in_bounds]
    if len(seam_points) == 0:
        raise RuntimeError('No seam points remain after raster bounds filter')

    if shared_edge is not None:
        orient, const, v0, v1 = shared_edge
        if orient == "v":
            sel_a = (np.abs(verts_a[:, 0] - const) <= snap_tol) & (verts_a[:, 1] >= v0) & (verts_a[:, 1] <= v1)
            sel_b = (np.abs(verts_b[:, 0] - const) <= snap_tol) & (verts_b[:, 1] >= v0) & (verts_b[:, 1] <= v1)
        else:
            sel_a = (np.abs(verts_a[:, 1] - const) <= snap_tol) & (verts_a[:, 0] >= v0) & (verts_a[:, 0] <= v1)
            sel_b = (np.abs(verts_b[:, 1] - const) <= snap_tol) & (verts_b[:, 0] >= v0) & (verts_b[:, 0] <= v1)
        if np.any(sel_a):
            seam_points = np.vstack([seam_points, verts_a[sel_a]])
        if np.any(sel_b):
            seam_points = np.vstack([seam_points, verts_b[sel_b]])

    point_cap = args.get('seam_point_cap', None)
    if point_cap is not None and len(seam_points) > point_cap:
        idx = np.linspace(0, len(seam_points) - 1, point_cap, dtype=int)
        seam_points = seam_points[idx]
        print(f'Seam point cap {point_cap} applied, kept {len(seam_points)} points')

    seam_points = dedupe_points(seam_points, snap_tol)

    points_file = args['out_prefix'] + '_points.txt'
    with open(points_file, 'w') as f:
        for v in seam_points:
            f.write(f'{v[0]} {v[1]}\n')

    seam_lloyd = args.get('seam_lloyd', None)
    if seam_lloyd is not None:
        args = dict(args)
        args['lloyd_itr'] = int(seam_lloyd)

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
    spatial_map = {}
    tris_out = []


    for tri in tris_a[~mask_a]:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, spatial_map, verts_a[idx], snap_tol))
        tris_out.append(idxs)

    for tri in tris_b[~mask_b]:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, spatial_map, verts_b[idx], snap_tol))
        tris_out.append(idxs)

    for tri in seam_tris:
        idxs = []
        for idx in tri:
            idxs.append(add_vertex(verts_out, index_map, spatial_map, seam_verts[idx], snap_tol))
        tris_out.append(idxs)

    verts_out = np.asarray(verts_out, dtype=float)
    tris_out = np.asarray(tris_out, dtype=int)

    merged_core_bbox = [
        min(bbox_a[0], bbox_b[0]),
        min(bbox_a[1], bbox_b[1]),
        max(bbox_a[2], bbox_b[2]),
        max(bbox_a[3], bbox_b[3])
    ]

    band_mask = compute_band_mask(verts_out, tris_out, merged_core_bbox, band_width)

    npz_path = args['out_prefix'] + '.npz'
    np.savez(npz_path, verts=verts_out, tris=tris_out, band_tri_mask=band_mask)

    meta_path = args['out_prefix'] + '.json'
    with open(meta_path, 'w') as f:
        json.dump({
            'tile_bbox': merged_core_bbox,
            'core_bbox': merged_core_bbox,
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
    if len(sys.argv) >= 2 and sys.argv[1] == '--single':
        if len(sys.argv) < 5:
            raise SystemExit('usage: MPI_merge_tiles.py --single tile_a_npz tile_b_npz out_prefix [meta_a] [meta_b]')
        tile_a_npz = sys.argv[2]
        tile_b_npz = sys.argv[3]
        out_prefix = sys.argv[4]
        tile_a_meta = sys.argv[5] if len(sys.argv) > 5 else tile_a_npz.replace('.npz', '.json')
        tile_b_meta = sys.argv[6] if len(sys.argv) > 6 else tile_b_npz.replace('.npz', '.json')

        args = {
            'tile_a': {'npz': tile_a_npz, 'meta': tile_a_meta},
            'tile_b': {'npz': tile_b_npz, 'meta': tile_b_meta},
            'out_prefix': out_prefix,
            'gdal_prefix': os.environ.get('GDAL_PREFIX', ''),
            'mesher_path': os.environ.get('MESHER_EXE', ''),
            'constraints': [],
            'parameter_files': {},
            'initial_conditions': {},
            'max_area': 1,
            'min_area': 1,
            'max_tolerance': 1,
            'errormetric': 'rmse',
            'lloyd_itr': 0,
            'use_weights': False,
            'topo_weight': 1.0,
            'weight_threshold': 0.0,
            'is_geographic': False,
            'dem_path': os.environ.get('MESHER_DEM', ''),
            'seam_point_spacing': None,
            'merge_snap_tol': None,
            'debug_seam': True
        }
        if not args['dem_path']:
            raise SystemExit('MESHER_DEM env var is required for --single')
        if not args['mesher_path']:
            raise SystemExit('MESHER_EXE env var is required for --single')
        merge_tiles(args)
    else:
        main(*sys.argv[1:])
