"""Replace artificial MPI boundaries with validated constrained mesh cavities.

Only whole triangles are retained. The patch is the exact set difference
between the intended domain and their coverage, never a centroid-clipped strip.
CGAL may refine the patch perimeter: if it changes a retained interface, discard
the candidate and enlarge the cavity. No tolerance-based vertex welding is used.
"""

import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from osgeo import ogr
from shapely import STRtree, from_wkb, polygons, union_all
from shapely.geometry import LineString, Polygon, box, mapping
from shapely.ops import polygonize

from mesher.mesh_topology import MeshTopologyError, build_neighbours, boundary_subdivisions

ogr.UseExceptions()


class PatchError(MeshTopologyError):
    """A proposed replacement cannot be accepted."""


class PatchQualityError(PatchError):
    def __init__(self, message, failures):
        self.failures = failures
        super().__init__(message)


class InterfaceChanged(PatchError):
    def __init__(self, faces):
        self.faces = np.unique(faces)
        super().__init__(f'Refinement changed a retained interface at {len(self.faces)} faces; cavity must grow')


def polygon_domain(geometry):
    """Keep all polygon components and holes, rejecting malformed domains."""
    if geometry.is_empty:
        return Polygon()
    if geometry.geom_type == 'GeometryCollection':
        geometry = union_all([g for g in geometry.geoms if g.geom_type in ('Polygon', 'MultiPolygon')])
    if geometry.geom_type not in ('Polygon', 'MultiPolygon') or not geometry.is_valid:
        raise PatchError('Domain must contain valid polygons; refusing implicit geometry repair')
    return geometry


def load_domain(path):
    ds = ogr.Open(str(path))
    if ds is None:
        raise PatchError(f'Cannot open domain {path}')
    pieces = [from_wkb(bytes(f.GetGeometryRef().ExportToWkb())) for f in ds.GetLayer()
              if f.GetGeometryRef() is not None]
    return polygon_domain(union_all(pieces))


def load_constraints(paths):
    lines = []
    for path in paths:
        ds = ogr.Open(str(path))
        if ds is None:
            raise PatchError(f'Cannot open constraint {path}')
        for feature in ds.GetLayer():
            geom = feature.GetGeometryRef()
            if geom is not None:
                lines.extend(line_parts(from_wkb(bytes(geom.ExportToWkb()))))
    return lines


def line_parts(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type in ('LineString', 'LinearRing'):
        return [LineString(geometry.coords)]
    if geometry.geom_type in ('Polygon', 'MultiPolygon'):
        return line_parts(geometry.boundary)
    if hasattr(geometry, 'geoms'):
        return [line for part in geometry.geoms for line in line_parts(part)]
    return []  # Point intersections do not define a segment.


def mesh_coverage(vertices, triangles):
    """Polygonize directed boundary edges without unioning all interior faces.

    A polygonized face is retained only if its CCW exterior follows a directed
    triangle boundary. This excludes holes and handles point-touching components.
    Input meshes must already be planar (CGAL output or a validated merge).
    """
    if not len(triangles):
        return Polygon()
    neighbours = build_neighbours(triangles, len(vertices), vertices)
    faces, slots = np.nonzero(neighbours < 0)
    edges = []
    directed = set()
    for face, slot in zip(faces, slots):
        tri = triangles[face]
        p, q, r = vertices[tri, :2]
        a, b = vertices[tri[(slot + 1) % 3], :2], vertices[tri[(slot + 2) % 3], :2]
        if (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]) < 0:
            a, b = b, a
        pair = (tuple(a), tuple(b))
        directed.add(pair)
        edges.append(LineString(pair))
    pieces = []
    for poly in polygonize(edges):
        coords = list(poly.exterior.coords)
        if not poly.exterior.is_ccw:
            coords.reverse()
        if (coords[0], coords[1]) in directed:
            pieces.append(poly)
    coverage = polygon_domain(union_all(pieces))
    pts = vertices[triangles, :2]
    u, v = pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0]
    area = .5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]).sum()
    if abs(area - coverage.area) > area_tolerance(coverage):
        raise PatchError('Retained mesh has inconsistent boundary coverage')
    return coverage


def area_tolerance(domain):
    if domain.is_empty:
        return 1e-10
    scale = max(1.0, *(abs(x) for x in domain.bounds))
    return max(1e-10, 64 * np.finfo(float).eps * scale * domain.length)


def validate_patch(vertices, triangles, domain, physical_constraints=()):
    """Check degeneracy, positive-area overlaps and exact intended coverage."""
    neighbours = build_neighbours(triangles, len(vertices), vertices)
    geometries = polygons(vertices[triangles, :2])
    coverage = union_all(geometries)
    epsilon = area_tolerance(domain)
    total_area = sum(g.area for g in geometries)
    if total_area - coverage.area > epsilon:
        raise PatchError('Replacement triangles overlap')
    if coverage.symmetric_difference(domain).area > epsilon:
        raise PatchError('Replacement triangles do not cover exactly the cavity domain')
    # Constraints must lie on triangle edges, not merely inside the covered area.
    if physical_constraints:
        edges = union_all([g.boundary for g in geometries])
        length_tolerance = max(1e-10, epsilon / max(1.0, domain.length))
        for line in physical_constraints:
            clipped = line.intersection(domain)
            if clipped.difference(edges).length > length_tolerance:
                raise PatchError('A physical constraint is missing from the replacement mesh')
    return neighbours


def _mesher_arguments(options):
    args = [options['mesher_path'], '--area', str(options['max_area']),
            '--min-area', str(options['min_area']), '--error-metric', options['errormetric'],
            '--raster', options['dem_path'], '--tolerance', str(options['max_tolerance'])]
    weighted = options.get('use_weights', False)
    if weighted:
        args += ['--weight', str(options['topo_weight']), '--weight-threshold', str(options['weight_threshold'])]
    continuous, categorical = [], []
    for files in (options.get('parameter_files', {}), options.get('initial_conditions', {})):
        for data in files.values():
            if 'tolerance' not in data:
                continue
            (categorical if data.get('method') == 'mode' else continuous).append(data)
    # The C++ reader consumes all continuous rasters before categorical rasters,
    # regardless of their command-line order. Weights must follow that order too.
    for records, raster_flag, tolerance_flag in ((continuous, '--raster', '--tolerance'),
                                               (categorical, '--category-raster', '--category-frac')):
        for data in records:
            filename = data['filename']
            if not isinstance(filename, str):
                filename = filename[0]
            args += [raster_flag, filename, tolerance_flag, str(data['tolerance'])]
            if weighted:
                args += ['--weight', str(data.get('weight', 1))]
    if options.get('mesher_debug', False):
        args += ['--debug', 'true']
    return args


def remesh_domain(domain, constraints, options, directory, fixed_edges=()):
    """Invoke the real C++ refinement criteria, including raster constraints."""
    domain = polygon_domain(domain)
    if domain.is_empty:
        return np.empty((0, 3)), np.empty((0, 3), dtype=np.int64)
    if options.get('is_geographic', False):
        raise PatchError('MPI remeshing requires a projected working CRS')
    # Node all linework together before passing it to CGAL. Never drop crossings.
    linework = [domain.boundary] + [line.intersection(domain) for line in constraints] + list(fixed_edges)
    noded = line_parts(union_all(linework))
    features = [{'type': 'Feature', 'properties': {}, 'geometry': mapping(line)} for line in noded]
    # Preserve failed inputs and the full backend log for production cases.
    # Stream output to disk: repeated GDAL errors must not fill the parent's RAM.
    work = Path(tempfile.mkdtemp(prefix='cavity-', dir=directory)).resolve()
    try:
        (work / 'domain.json').write_text(json.dumps(mapping(domain)))
        (work / 'constraints.json').write_text(json.dumps({'type': 'FeatureCollection', 'features': features}))
        output = work / 'mesh.json'
        command = _mesher_arguments(options) + ['--domain-file', str(work / 'domain.json'),
                   '--interior-plgs-file', str(work / 'constraints.json'), '--patch-output', str(output),
                   '--patch-max-vertices', str(options.get('remesh_max_vertices', 500000))]
        (work / 'invocation.json').write_text(json.dumps({'command': command, 'cwd': os.getcwd(),
            'tile': options.get('tile_prefix'), 'merge': options.get('out_prefix'),
            'domain_bounds': list(domain.bounds), 'domain_area': domain.area}, indent=2))
        with (work / 'backend.log').open('w') as log:
            run = subprocess.run(command, text=True, stdout=log, stderr=subprocess.STDOUT)
        if run.returncode:
            with (work / 'backend.log').open('rb') as log:
                log.seek(0, os.SEEK_END)
                log.seek(max(0, log.tell() - 3000))
                tail = log.read().decode(errors='replace')
            message = f'C++ patch refinement failed (exit {run.returncode}): ' + tail
            quality_path = Path(str(output) + '.quality.json')
            if quality_path.exists():
                raise PatchQualityError(message, json.loads(quality_path.read_text())['faces'])
            raise PatchError(message)
        data = json.loads(output.read_text())
        vertices = np.asarray(data['vertices'], dtype=float)
        triangles = np.asarray(data['triangles'], dtype=np.int64)
        validate_patch(vertices, triangles, domain, constraints)
    except Exception as exc:
        context = options.get('tile_prefix') or options.get('out_prefix') or 'cavity'
        message = f'{context}: {exc}\nReproduction files and full backend log: {work}'
        if isinstance(exc, PatchError):
            exc.args = (message,)
            raise
        raise PatchError(message) from exc
    shutil.rmtree(work)
    return vertices, triangles


def _shared_interface(a, b):
    """Return the exact coordinate of a shared rectangular core interface."""
    for axis in (0, 1):
        for first, second in ((a, b), (b, a)):
            if first[axis + 2] == second[axis]:
                other = 1 - axis
                if min(first[other + 2], second[other + 2]) > max(first[other], second[other]):
                    return axis, first[axis + 2]
    raise PatchError('Tile core boxes must share an exact edge')


def _retain_mask(vertices, triangles, bounds, axis, seam, width):
    pts = vertices[triangles, :2]
    inside = ((pts >= np.asarray(bounds[:2])) & (pts <= np.asarray(bounds[2:]))).all(axis=(1, 2))
    # All three vertices must be on the retained side; no centroid clipping.
    away = (np.abs(pts[:, :, axis] - seam) >= width).all(axis=1)
    return inside & away


def _quality_interface_faces(fixed_edges, owners, failures, domain):
    """Locate retained faces touched by bad patch triangles, without welding.

    A small coordinate-roundoff allowance handles refined points on a sloping
    edge. It is used only to select whole faces for removal, never to accept a
    candidate, move vertices or alter the physical domain.
    """
    if not fixed_edges or not failures:
        return np.empty(0, dtype=np.int64)
    tree = STRtree(fixed_edges)
    scale = max(1.0, *(abs(x) for x in domain.bounds))
    epsilon = max(1e-10, 64 * np.finfo(float).eps * scale)
    affected = set()
    for failure in failures:
        triangle = Polygon(failure['vertices'])
        x0, y0, x1, y1 = triangle.bounds
        candidates = tree.query(box(x0-epsilon, y0-epsilon, x1+epsilon, y1+epsilon))
        for edge in candidates:
            if fixed_edges[edge].distance(triangle) <= epsilon:
                affected.add(owners[edge])
    return np.asarray(sorted(affected), dtype=np.int64)


def _stitch(retained_v, retained_t, patch_v, patch_t, domain):
    # Only exact XY equality to a retained vertex authorizes reusing its ID.
    # The C++ writer uses max_digits10 so interface doubles round-trip unchanged.
    used = np.unique(retained_t)
    vertices = retained_v[used].tolist()
    mapping_old = np.full(len(retained_v), -1, dtype=np.int64)
    mapping_old[used] = np.arange(len(used))
    kept = mapping_old[retained_t]
    ids = {tuple(v[:2]): i for i, v in enumerate(vertices)}
    if len(ids) != len(vertices):
        raise PatchError('Distinct retained vertices have coincident coordinates')
    remap = np.empty(len(patch_v), dtype=np.int64)
    for i, point in enumerate(patch_v):
        key = tuple(point[:2])
        if key not in ids:
            ids[key] = len(vertices)
            vertices.append(point.tolist())
        remap[i] = ids[key]
    result_v = np.asarray(vertices, dtype=float)
    result_t = np.vstack((kept, remap[patch_t]))
    neighbours = build_neighbours(result_t, len(result_v), result_v)
    # Every retained edge which faces the cavity must have exactly one patch
    # neighbour. Boundary points that CGAL inserted therefore trigger a retry.
    if len(kept):
        old_neighbours = build_neighbours(kept, len(result_v), result_v)
        faces, slots = np.nonzero(old_neighbours < 0)
        changed = []
        for face, slot in zip(faces, slots):
            edge = result_v[kept[face, [(slot + 1) % 3, (slot + 2) % 3]], :2]
            line = LineString(edge)
            if not domain.boundary.covers(line) and neighbours[face, slot] < len(kept):
                changed.append(face)
        if changed:
            raise InterfaceChanged(changed)
    coincident, hanging = boundary_subdivisions(result_v, result_t, neighbours, 1e-8)
    if coincident or next(hanging, None) is not None:
        raise PatchError('Joined mesh has coincident or hanging boundary nodes')
    return result_v, result_t


def merge_cavity(vertices_a, triangles_a, bounds_a, vertices_b, triangles_b, bounds_b,
                 domain, constraints, options, directory, mesher=remesh_domain):
    """Retry a local cavity transaction; never modify either input tile."""
    axis, seam = _shared_interface(bounds_a, bounds_b)
    bounds = [min(bounds_a[0], bounds_b[0]), min(bounds_a[1], bounds_b[1]),
              max(bounds_a[2], bounds_b[2]), max(bounds_a[3], bounds_b[3])]
    domain = polygon_domain(domain.intersection(box(*bounds)))
    if domain.is_empty:
        return np.empty((0, 3)), np.empty((0, 3), dtype=np.int64), {'attempts': 0}
    supplied_width = options.get('remesh_band_width')
    width = float(2 * math.sqrt(options['max_area']) if supplied_width is None else supplied_width)
    attempts = int(options.get('remesh_max_attempts', 12))
    if not np.isfinite(width) or width <= 0 or attempts < 1:
        raise PatchError('Remesh band width and attempt count must be positive')
    joined_v = np.vstack((vertices_a, vertices_b))
    joined_t = np.vstack((triangles_a, triangles_b + len(vertices_a)))
    keep = np.concatenate((_retain_mask(vertices_a, triangles_a, bounds_a, axis, seam, width),
                           _retain_mask(vertices_b, triangles_b, bounds_b, axis, seam, width)))
    # Tile-local neighbours are sufficient for growing the cavity into either
    # retained core. Overlapping halo vertices still have separate IDs here.
    original_neighbours = build_neighbours(joined_t, len(joined_v), joined_v) if len(joined_t) else None
    last_error = None
    for attempt in range(attempts):
        retained_ids = np.flatnonzero(keep)
        retained = joined_t[keep]
        coverage = mesh_coverage(joined_v, retained)
        if coverage.difference(domain).area > area_tolerance(domain):
            raise PatchError('Retained triangles extend outside the intended domain')
        cavity = polygon_domain(domain.difference(coverage))
        print(f'Remesh attempt {attempt + 1}: width={width:g}, retained={len(retained)}, cavity_area={cavity.area:g}', flush=True)
        try:
            # Preserve every retained boundary vertex, including collinear ones
            # which a polygon boolean operation might otherwise discard.
            fixed, fixed_owners = [], []
            if len(retained):
                neighbours = build_neighbours(retained, len(joined_v), joined_v)
                faces, slots = np.nonzero(neighbours < 0)
                for f, s in zip(faces, slots):
                    line = LineString(joined_v[retained[f, [(s + 1) % 3, (s + 2) % 3]], :2])
                    if cavity.boundary.covers(line) and not domain.boundary.covers(line):
                        fixed.append(line)
                        fixed_owners.append(f)
            patch_v, patch_t = mesher(cavity, constraints, options, directory, fixed)
            if len(patch_t):
                validate_patch(patch_v, patch_t, cavity, constraints)
            output_v, output_t = _stitch(joined_v, retained, patch_v, patch_t, domain)
            # Coverage of retained triangles and patch was checked separately;
            # set difference makes their interiors disjoint by construction.
            return output_v, output_t, {'attempts': attempt + 1, 'band_width': width,
                                       'retained_triangles': len(retained), 'patch_triangles': len(patch_t)}
        except (PatchError, MeshTopologyError) as exc:
            last_error = exc
            print(f'Remesh candidate rejected: {exc}', flush=True)
            local_faces = np.empty(0, dtype=np.int64)
            if isinstance(exc, InterfaceChanged):
                local_faces = exc.faces
            elif isinstance(exc, PatchQualityError):
                local_faces = _quality_interface_faces(fixed, fixed_owners, exc.failures, domain)
            if len(local_faces):
                # Both split interfaces and unresolved quality failures can be
                # caused by retained protrusions. Remove the touched faces and
                # one neighbour ring, leaving unrelated seam regions alone.
                affected = retained_ids[local_faces]
                adjacent = original_neighbours[affected].ravel()
                keep[affected] = False
                keep[adjacent[adjacent >= 0]] = False
                print(f'Growing cavity locally at {len(affected)} retained faces plus neighbours', flush=True)
                continue
            if not len(retained):
                break
            if 'Patch vertex budget exceeded' in str(exc):
                break
            width *= 2
            keep &= np.concatenate((_retain_mask(vertices_a, triangles_a, bounds_a, axis, seam, width),
                                    _retain_mask(vertices_b, triangles_b, bounds_b, axis, seam, width)))
    raise PatchError(f'No valid cavity replacement after {attempt + 1} attempts: {last_error}')
