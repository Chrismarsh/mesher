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
import subprocess
import numpy as np
from mpi4py import MPI
from osgeo import ogr, gdal
from mesher.mesher_utls.ogr_utils import load_polygon_geom, normalize_polygon, \
    linestring_features_from_geom, extract_polygons

gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()


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


def write_polygon_geojson(path, coords):
    if len(coords) < 4:
        raise RuntimeError('Cannot write seam GeoJSON with fewer than 4 points')
    if coords[0] != coords[-1]:
        coords = coords + [coords[0]]
    gj = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {},
            "geometry": {"type": "Polygon", "coordinates": [coords]}
        }]
    }
    with open(path, 'w') as f:
        json.dump(gj, f)


def write_points_geojson(path, points):
    features = []
    for x, y in points:
        features.append({
            "type": "Feature",
            "properties": {},
            "geometry": {"type": "Point", "coordinates": [float(x), float(y)]}
        })
    gj = {"type": "FeatureCollection", "features": features}
    with open(path, 'w') as f:
        json.dump(gj, f)


def densify_ring_coords(coords, max_len):
    if len(coords) < 2 or max_len <= 0:
        return coords

    densified = [coords[0]]
    for i in range(1, len(coords)):
        x0, y0 = densified[-1]
        x1, y1 = coords[i]
        dx = x1 - x0
        dy = y1 - y0
        seg_len = math.hypot(dx, dy)
        if seg_len <= max_len or seg_len == 0:
            densified.append([x1, y1])
            continue
        steps = int(math.ceil(seg_len / max_len))
        for k in range(1, steps + 1):
            t = float(k) / float(steps)
            densified.append([x0 + t * dx, y0 + t * dy])

    if densified[0] != densified[-1]:
        densified.append(densified[0])
    return densified


def median_edge_length(verts, tris, mask):
    if len(tris) == 0 or not np.any(mask):
        return None
    lengths = []
    for tri in tris[mask]:
        v0 = verts[tri[0]]
        v1 = verts[tri[1]]
        v2 = verts[tri[2]]
        lengths.append(math.hypot(v1[0] - v0[0], v1[1] - v0[1]))
        lengths.append(math.hypot(v2[0] - v1[0], v2[1] - v1[1]))
        lengths.append(math.hypot(v0[0] - v2[0], v0[1] - v2[1]))
    if not lengths:
        return None
    return float(np.median(lengths))


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
    if args.get('mesher_debug', False):
        execstr += ' --debug true'

    if args.get('skip_angle_below_min_area', False):
        execstr += ' --skip-angle-below-min-area true'

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


def triangle_intersects_geometry(verts, tris, geom):
    mask = np.zeros(len(tris), dtype=bool)
    if len(tris) == 0:
        return mask

    poly_env = geom.GetEnvelope()
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

        if poly_tri.Intersects(geom):
            mask[i] = True

    return mask


def triangle_intersects_polygon(verts, tris, poly):
    return triangle_intersects_geometry(verts, tris, poly)


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

    if meta_a.get('empty_tile') or meta_b.get('empty_tile'):
        out_prefix = args['out_prefix']
        out_npz = out_prefix + '.npz'
        out_meta = out_prefix + '.json'
        if meta_a.get('empty_tile') and meta_b.get('empty_tile'):
            verts = np.zeros((0, 3), dtype=float)
            tris = np.zeros((0, 3), dtype=int)
            band_mask = np.zeros((0,), dtype=bool)
            np.savez(out_npz, verts=verts, tris=tris, band_tri_mask=band_mask)
            with open(out_meta, 'w') as f:
                json.dump({
                    'tile_bbox': meta_a.get('tile_bbox', meta_b.get('tile_bbox')),
                    'core_bbox': meta_a.get('core_bbox', meta_b.get('core_bbox')),
                    'band_width': meta_a.get('band_width', meta_b.get('band_width')),
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

    verts_a = data_a['verts']
    tris_a = data_a['tris']
    band_a = data_a['band_tri_mask']
    verts_b = data_b['verts']
    tris_b = data_b['tris']
    band_b = data_b['band_tri_mask']

    log_step('Computing seam bbox/strip')
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

    log_step('Clipping seam to raster bounds')
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

    seam_strip = normalize_polygon(seam_strip.Intersection(raster_poly))
    if seam_strip is None:
        raise RuntimeError('Seam strip clipped outside raster bounds')

    log_step('Selecting seam/ownership masks')
    # remove any triangles from either tile that are in the seam strip or outside ownership
    # ownership is centroid-based; only enforce outside-core pruning within the seam zone
    inside_a = centroid_mask_in_bbox(verts_a, tris_a, bbox_a)
    inside_b = centroid_mask_in_bbox(verts_b, tris_b, bbox_b)

    # initial seam selection based on strip around shared edge
    seam_select_a = triangle_intersects_polygon(verts_a, tris_a, seam_strip)
    seam_select_b = triangle_intersects_polygon(verts_b, tris_b, seam_strip)

    outside_a = (~inside_a) & centroid_mask_in_bbox(verts_a, tris_a, seam_bbox)
    outside_b = (~inside_b) & centroid_mask_in_bbox(verts_b, tris_b, seam_bbox)

    # removal mask includes seam strip and non-owned triangles near the seam
    mask_a = seam_select_a | outside_a
    mask_b = seam_select_b | outside_b
    if not np.any(mask_a):
        mask_a = band_a
    if not np.any(mask_b):
        mask_b = band_b

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

    snap_tol = args.get('merge_snap_tol', None)
    seam_spacing = args.get('seam_point_spacing', None)
    if snap_tol is None:
        if seam_spacing is None:
            snap_tol = 1e-6
        else:
            snap_tol = max(1e-6, seam_spacing * 0.01)

    # Build an irregular seam polygon from the triangles we will remove.
    log_step('Building seam polygon from removed triangles')
    if np.any(mask_a) or np.any(mask_b):
        log_step('Union triangles from tile A')
        seam_geom = triangles_to_union_polygon(verts_a, tris_a, mask_a)
        log_step('Union triangles from tile B')
        seam_geom = seam_geom.Union(triangles_to_union_polygon(verts_b, tris_b, mask_b))
        log_step('Extracting seam polygon')
        seam_poly = extract_polygons(seam_geom)
        seam_poly = seam_poly[0] if seam_poly else None
    else:
        seam_poly = seam_strip

    outer_polygon_shp = args.get('outer_polygon_shp', None)
    outer_poly = None
    if outer_polygon_shp:
        outer_poly = load_polygon_geom(outer_polygon_shp)

    if seam_spacing is not None:
        log_step('Buffering seam for removal expansion')
        buf_dist = min(float(seam_spacing) * 0.25, band_width * 0.5)
        if buf_dist > 0:
            buffered = seam_poly.Buffer(buf_dist)
            buffered = normalize_polygon(buffered.Intersection(raster_poly))
            if buffered is not None and outer_poly is not None:
                buffered = normalize_polygon(buffered.Intersection(outer_poly))
            if buffered is None:
                raise RuntimeError('Buffered seam polygon invalid after clip')
            # Expand removal to match the buffered seam area to avoid overlaps.
            mask_a = mask_a | triangle_intersects_polygon(verts_a, tris_a, buffered)
            mask_b = mask_b | triangle_intersects_polygon(verts_b, tris_b, buffered)
            seam_geom = triangles_to_union_polygon(verts_a, tris_a, mask_a)
            seam_geom = seam_geom.Union(triangles_to_union_polygon(verts_b, tris_b, mask_b))
            seam_poly = normalize_polygon(seam_geom)
            if seam_poly is None:
                raise RuntimeError('Seam polygon invalid after buffer expansion')

    log_step('Clipping seam to raster bounds (final)')
    seam_poly = normalize_polygon(seam_poly.Intersection(raster_poly))
    if seam_poly is None:
        raise RuntimeError('Seam polygon invalid after raster bounds clipping')

    if outer_poly is not None:
        log_step('Clipping seam to outer polygon')
        seam_poly = normalize_polygon(seam_poly.Intersection(outer_poly))
        if seam_poly is None:
            raise RuntimeError('Seam polygon invalid after outer polygon clipping')

    # Guard against skinny seam polygons which can cause the mesher to hang.
    seam_env = seam_poly.GetEnvelope()
    seam_w = seam_env[1] - seam_env[0]
    seam_h = seam_env[3] - seam_env[2]
    if seam_w <= 0 or seam_h <= 0:
        raise RuntimeError('Seam polygon has non-positive bounds')
    seam_aspect = max(seam_w / seam_h, seam_h / seam_w)
    aspect_limit = float(args.get('seam_aspect_limit', 10.0))
    min_edge_limit = 1.25 * band_width
    if seam_aspect > aspect_limit and min(seam_w, seam_h) <= min_edge_limit:
        raise RuntimeError(
            f'Seam polygon aspect ratio {seam_aspect:.2f} exceeds limit {aspect_limit:.2f}. '
            f'Seam bounds={seam_env}.')

    log_step('Aligning removal masks to final seam polygon')
    # Align removal to the final seam polygon to avoid holes.
    seam_select_a = triangle_intersects_polygon(verts_a, tris_a, seam_poly)
    seam_select_b = triangle_intersects_polygon(verts_b, tris_b, seam_poly)
    mask_a = mask_a & seam_select_a
    mask_b = mask_b & seam_select_b
    if not np.any(mask_a):
        mask_a = seam_select_a
    if not np.any(mask_b):
        mask_b = seam_select_b

    if args.get('seam_constraints', True):
        log_step('Clipping constraints to seam')
        for cpath in args['constraints']:
            ds = ogr.Open(cpath)
            if ds is None:
                raise RuntimeError(f'Unable to open constraint {cpath}')
            layer = ds.GetLayer(0)
            for feat in layer:
                geom = feat.GetGeometryRef()
                if geom is None:
                    continue
                clipped = geom.Intersection(seam_poly)
                if clipped is None:
                    continue
                interior_PLGS['features'].extend(linestring_features_from_geom(clipped))
            ds = None

    interior_plgs_file = args['out_prefix'] + '_interior_PLGS.geojson'
    log_step('Writing interior PLGS geojson')
    with open(interior_plgs_file, 'w') as fp:
        json.dump(interior_PLGS, fp)

    log_step('Building seam poly points')
    coords = polygon_exterior_coords(seam_poly)
    seam_simplify_tol = args.get('seam_simplify_tol', None)
    if seam_simplify_tol is not None:
        coords = clean_ring_coords(coords, float(seam_simplify_tol))
        if len(coords) < 4:
            raise RuntimeError('Seam polygon simplified to too few points')
    if seam_spacing is not None:
        log_step('Densifying seam boundary')
        seam_boundary = seam_poly.GetBoundary()
        boundary_a = triangle_intersects_geometry(verts_a, tris_a, seam_boundary)
        boundary_b = triangle_intersects_geometry(verts_b, tris_b, seam_boundary)
        edge_a = median_edge_length(verts_a, tris_a, boundary_a)
        edge_b = median_edge_length(verts_b, tris_b, boundary_b)
        edge_lengths = [v for v in (edge_a, edge_b) if v is not None]
        if edge_lengths:
            boundary_spacing = float(np.median(edge_lengths))
        else:
            boundary_spacing = min(float(seam_spacing) * 0.5, band_width)
        coords = densify_ring_coords(coords, boundary_spacing)
    poly_file = args['out_prefix'] + '_seam.poly'
    write_poly_from_coords(poly_file, coords)
    if args.get('dump_poly_files', False) or args.get('dump_poly_only', False):
        seam_geojson = args['out_prefix'] + '_seam.geojson'
        write_polygon_geojson(seam_geojson, coords)

    seam_points = []
    for x, y in coords:
        seam_points.append([x, y, 0.0])
    seam_points = np.asarray(seam_points, dtype=float)

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

    log_step('Dedupe seam points')
    seam_points = dedupe_points(seam_points, snap_tol)

    points_file = args['out_prefix'] + '_points.txt'
    with open(points_file, 'w') as f:
        for v in seam_points:
            f.write(f'{v[0]} {v[1]}\n')
    if args.get('dump_poly_files', False) or args.get('dump_poly_only', False):
        points_geojson = args['out_prefix'] + '_points.geojson'
        write_points_geojson(points_geojson, seam_points[:, :2])

    if args.get('dump_poly_only', False):
        print(f'Dumping seam poly only (no mesher run): {poly_file}')
        return
    if args.get('dump_poly_files', False):
        print(f'Dumping seam poly files (mesher will still run): {poly_file}')

    seam_lloyd = args.get('seam_lloyd', None)
    if seam_lloyd is not None:
        args = dict(args)
        args['lloyd_itr'] = int(seam_lloyd)

    execstr = build_mesher_exec_str(args, poly_file, interior_plgs_file, points_file)
    print(f'Running mesher: {execstr}', flush=True)
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
            raise SystemExit('usage: MPI_merge_tiles.py --single tile_a_npz tile_b_npz out_prefix [meta_a] [meta_b] [--dump-only]')
        tile_a_npz = sys.argv[2]
        tile_b_npz = sys.argv[3]
        out_prefix = sys.argv[4]
        dump_only = '--dump-only' in sys.argv
        args_in = [a for a in sys.argv if a != '--dump-only']
        tile_a_meta = args_in[5] if len(args_in) > 5 else tile_a_npz.replace('.npz', '.json')
        tile_b_meta = args_in[6] if len(args_in) > 6 else tile_b_npz.replace('.npz', '.json')

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
            'dump_poly_only': dump_only,
            'dump_poly_files': dump_only
        }
        if not args['dem_path']:
            raise SystemExit('MESHER_DEM env var is required for --single')
        if not args['mesher_path']:
            raise SystemExit('MESHER_EXE env var is required for --single')
        merge_tiles(args)
    else:
        main(*sys.argv[1:])
