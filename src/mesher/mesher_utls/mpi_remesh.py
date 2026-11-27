"""Replace artificial MPI boundaries with validated constrained mesh cavities.

Only whole triangles are retained. The patch is the exact set difference
between the intended domain and their coverage, never a centroid-clipped strip.
CGAL may refine patch interfaces. Subdivide their retained owners first, then
repair local two-sided regions while preserving validated triangles elsewhere.
No tolerance-based vertex welding is used.
"""

import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import numpy as np
from osgeo import ogr
from shapely import (STRtree, box as boxes, covers, distance, from_wkb, get_coordinates,
                     linestrings, points, polygons, prepare, union_all)
from shapely.geometry import LineString, Point, Polygon, box, mapping
from shapely.geometry.polygon import orient
from shapely.ops import polygonize

from mesher.mesh_topology import MeshTopologyError, build_neighbours, boundary_subdivisions, boundary_node_conflicts
from mesher.mesher_utls.mpi_progress import debug_log, report_phase

ogr.UseExceptions()


class PatchError(MeshTopologyError):
    """A proposed replacement cannot be accepted."""


class BackendTerminated(PatchError):
    """A killed backend cannot be recovered by increasing the cavity."""


class PatchResult(tuple):
    """Validated immutable arrays plus explicit CGAL interface lineage."""
    # Keep two-value unpacking for mesh callers while carrying interface lineage and statistics
    # into repair. Read-only geometry can safely be reused by a retry.
    def __new__(cls, vertices, triangles, interfaces, stats, constrained_edges=None):
        vertices.setflags(write=False)
        triangles.setflags(write=False)
        result = super().__new__(cls, (vertices, triangles))
        result.interfaces = interfaces
        result.stats = stats
        result.constrained_edges = constrained_edges
        return result


def _run_backend(command, work, options):
    """Keep full output on disk and report progress without accumulating logs."""
    started = time.perf_counter()
    context = Path(options.get('out_prefix') or options.get('tile_prefix') or work).name
    with (work / 'backend.log').open('w') as log:
        with subprocess.Popen(command, text=True, stdout=log, stderr=subprocess.STDOUT) as process:
            while True:
                try:
                    code = process.wait(timeout=10)
                    break
                except subprocess.TimeoutExpired:
                    if not options.get('mesher_debug', False):
                        continue
                    with (work / 'backend.log').open('rb') as progress:
                        progress.seek(0, os.SEEK_END)
                        progress.seek(max(0, progress.tell() - 512))
                        lines = progress.read().decode(errors='replace').splitlines()
                    print(f'Backend active [{context}]: elapsed={time.perf_counter()-started:.1f}s; '
                          + (lines[-1] if lines else 'initializing'), flush=True)
    if code < 0:
        error = BackendTerminated(f'C++ backend terminated by signal {-code}; refusing cavity retry')
        error.stats = {'backend_seconds': time.perf_counter()-started}
        raise error
    return code, time.perf_counter() - started


class PatchQualityError(PatchError):
    # Keep the backend's per-face geometry and failure reasons so recovery can find the retained
    # owners near the failed faces.
    def __init__(self, message, failures):
        self.failures = failures
        super().__init__(message)


class PhysicalConstraintAngle(PatchQualityError):
    """An in-domain sector bounded by two immutable physical features."""


class PatchBoundaryRoundoff(PatchError):
    """Carry the unresolved geometry to the layer that owns its boundary faces."""
    # Pass the suspect chords and any known owners to the layer able to release them. Remember
    # whether the diagnosis came before native refinement.
    def __init__(self, message, chords, owners=None, *, preflight_only=False):
        self.preflight_only = preflight_only
        self.chords = chords
        self.owners = owners
        super().__init__(message)


class PatchEdgeError(PatchError):
    # Preserve the failed segment's endpoints so the caller can locate its owners in this
    # attempt's mesh.
    def __init__(self, message, edge):
        self.edge = edge
        super().__init__(message)


class PatchPrecisionError(PatchError):
    # Several near misses may identify the same retained interface. Carry each interface ID once
    # so the caller can release the corresponding faces.
    def __init__(self, interface_ids, epsilon):
        self.interface_ids = np.unique(interface_ids)
        super().__init__(f'Retained interfaces have a nonzero clearance below coordinate precision '
                         f'({epsilon:.6g}); cavity must grow at {len(self.interface_ids)} edges')


class InterfaceChanged(PatchError):
    """Face IDs index the retained_t argument of the failing _stitch call."""
    # Collect unique owners in retained-face numbering. The caller maps them back to
    # working-mesh rows before expanding the cavity.
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


# Combine the layer's polygon features into one validated footprint while keeping disconnected
# components and interior holes.
def load_domain(path):
    ds = ogr.Open(str(path))
    if ds is None:
        raise PatchError(f'Cannot open domain {path}')
    pieces = [from_wkb(bytes(f.GetGeometryRef().ExportToWkb())) for f in ds.GetLayer()
              if f.GetGeometryRef() is not None]
    return polygon_domain(union_all(pieces))


def _acute_corners(points):
    """Convex corners below CGAL's squared-sine bound, with interior on the left."""
    a, b = np.roll(points, 1, axis=0), np.roll(points, -1, axis=0)
    u, v = a-points, b-points
    lengths = np.linalg.norm(u, axis=1), np.linalg.norm(v, axis=1)
    cross = u[:, 0]*v[:, 1] - u[:, 1]*v[:, 0]
    with np.errstate(divide='ignore', invalid='ignore'):
        acute = (cross < 0) & (np.sum(u*v, axis=1) > 0) & ((cross / lengths[0] / lengths[1])**2 < .125)
    return [(points[i], a[i], b[i], min(lengths[0][i], lengths[1][i]))
            for i in np.flatnonzero(acute)]


def _acute_clip_corners(clipped, domain):
    """Find convex corners below the backend shape bound with an artificial edge."""
    candidates = []
    parts = [clipped] if clipped.geom_type == 'Polygon' else clipped.geoms
    for part in parts:
        part = orient(part)  # Domain interior lies to the left of every ring.
        for ring in [part.exterior, *part.interiors]:
            for p, a, b, length in _acute_corners(np.asarray(ring.coords[:-1])):
                candidates.append((p, a, b, length/4))
    if not candidates:
        return []
    epsilon = _coordinate_tolerance(domain)
    x0, y0, x1, y1 = clipped.bounds
    # Buffer only nearby physical linework, not the complete large domain.
    physical = domain.boundary.intersection(box(x0-epsilon, y0-epsilon, x1+epsilon, y1+epsilon)).buffer(epsilon)
    prepare(physical)
    return [(p, radius) for p, a, b, radius in candidates
            if not (physical.covers(LineString([a, p])) and physical.covers(LineString([p, b])))]


def _clipped_domain(domain, bounds, padding, corner_radius, context=None, mesher_debug=False):
    """Cover the entire owned core, extending only artificial acute clipping tips.

    A disk centred on a cut/physical-edge intersection makes the new cut meet
    that edge roughly perpendicularly. Only the clipping window grows; the
    physical domain, holes and every positive-area owned piece remain intact.
    The final group contains the whole domain and needs no artificial cuts.
    """
    core = box(*bounds)
    if domain.intersection(core).area == 0:
        return Polygon()
    window = box(bounds[0]-padding, bounds[1]-padding, bounds[2]+padding, bounds[3]+padding)

    # A halo can reach a disconnected island owned entirely by another tile. Include a clipped
    # component only if it has positive area inside this core.
    def owned_parts(window):
        clipped = polygon_domain(domain.intersection(window))
        parts = [clipped] if clipped.geom_type == 'Polygon' else clipped.geoms
        # Point/line-only contact owns nothing. Never apply a small-area cutoff.
        return polygon_domain(union_all([p for p in parts if p.intersection(core).area > 0]))

    clipped = owned_parts(window)
    initial_area, expanded = clipped.area, 0
    # Bounded, local geometry work; unresolved physical/constraint angles still
    # go through strict backend validation and preserve its reproduction files.
    for _ in range(4):
        corners = _acute_clip_corners(clipped, domain)
        disks = [Point(p).buffer(min(corner_radius, radius), quad_segs=4)
                 for p, radius in corners if min(corner_radius, radius) > _coordinate_tolerance(domain)]
        if not disks:
            break
        window = union_all([window, *disks])
        candidate = owned_parts(window)
        if candidate.equals(clipped):
            break
        expanded += len(disks)
        clipped = candidate
    if expanded and mesher_debug:
        label = f' [{Path(context).name}]' if context else ''
        print(f'Extended artificial clipping corners{label}: corners={expanded}, '
              f'added_halo_area={clipped.area-initial_area:.6g}', flush=True)
    return clipped


def tile_domain(domain, bounds, halo_width, context=None, mesher_debug=False):
    """Mesh owned halo components without introducing unavoidable cut angles."""
    return _clipped_domain(domain, bounds, halo_width, halo_width/2, context, mesher_debug)


# Read all configured feature layers into line segments that can be clipped to the domain of
# each tile or patch.
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


# Clipping can produce polygons, multipart lines and isolated points. Keep the linear pieces,
# including polygon rings, that can constrain a triangulation.
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


def mesh_coverage(vertices, triangles, neighbours=None):
    """Polygonize directed boundary edges without unioning all interior faces.

    A polygonized face is retained only if its CCW exterior follows a directed
    triangle boundary. This excludes holes and handles point-touching components.
    Input meshes must already be planar (CGAL output or a validated merge).
    """
    if not len(triangles):
        return Polygon()
    if neighbours is None:
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


# A roundoff-width displacement along a long perimeter can have measurable area. Scale coverage
# comparisons by perimeter and coordinate magnitude.
def area_tolerance(domain):
    if domain.is_empty:
        return 1e-10
    # At large projected eastings or northings, coordinate roundoff along the
    # perimeter can look like a thin strip of missing area. Perimeter times
    # coordinate error gives the area scale for this footprint check.
    scale = max(1.0, *(abs(x) for x in domain.bounds))
    return max(1e-10, 64 * np.finfo(float).eps * scale * domain.length)


# Coordinate precision depends on the magnitude of the projected values. Use this length scale
# for incidence checks that are sensitive to floating-point roundoff.
def _coordinate_tolerance(domain):
    scale = max(1.0, *(abs(x) for x in domain.bounds))
    return max(1e-10, 64 * np.finfo(float).eps * scale)


def _roundoff_chords(cavity, *, acute_only=False):
    """Find standalone rings whose transverse displacement fits one coordinate ULP."""
    components = [cavity] if cavity.geom_type == 'Polygon' else cavity.geoms
    chords = []
    for part in components:
        if part.is_empty or part.interiors:
            continue
        ulp = float(np.spacing(max(map(abs, part.bounds))))
        if part.area > part.length * ulp:
            continue
        xy = np.asarray(part.exterior.coords[:-1])
        axis = np.ptp(xy, axis=0).argmax()
        start, end = xy[xy[:, axis].argmin()], xy[xy[:, axis].argmax()]
        dx, dy = end - start
        offsets = xy - start
        if np.max(np.abs(dx * offsets[:, 1] - dy * offsets[:, 0])) <= math.hypot(dx, dy) * ulp:
            # With no mesh owner, refuse only a pinned corner that cannot
            # meet the shape bound. A tiny well-shaped physical island can.
            if not acute_only or _acute_corners(np.asarray(orient(part).exterior.coords[:-1])):
                chords.append(LineString([start, end]))
    return chords


def _unbisectable_edges(lines):
    """Return distinct endpoints whose computed midpoint rounds to an endpoint."""
    xy, line = get_coordinates(lines, return_index=True)
    a, b = xy[:-1], xy[1:]
    midpoint = a + (b-a)*.5
    collapsed = ((midpoint == a).all(axis=1) | (midpoint == b).all(axis=1))
    select = (line[:-1] == line[1:]) & (a != b).any(axis=1) & collapsed
    return list(linestrings(np.stack((a[select], b[select]), axis=1)))


def _roundoff_owners(chords, edges):
    """Find boundary edges following a chord over a shared longitudinal interval."""
    # An overlay-created endpoint can round off the original retained edge.
    # Use this allowance solely to choose whole faces to remesh.
    # A midpoint can round onto an endpoint at a junction of independently
    # cut edges. Select the incident rim faces at that point. For long
    # slivers, continue to require a shared interval along the chord.
    tolerances = [2 * np.spacing(max(map(abs, chord.bounds))) for chord in chords]
    point_contacts = [chord.length <= tol for chord, tol in zip(chords, tolerances)]
    envelopes = []
    for chord, tol in zip(chords, tolerances):
        x0, y0, x1, y1 = chord.bounds
        envelopes.append(box(x0-tol, y0-tol, x1+tol, y1+tol))
    tree = STRtree(envelopes)
    owners = set()
    for face, coordinates in edges:
        if face in owners:
            continue
        edge = LineString(coordinates)
        for i in tree.query(edge):
            chord, tol = chords[i], tolerances[i]
            if point_contacts[i]:
                if edge.distance(chord) <= tol:
                    owners.add(face)
                    break
                continue
            lo, hi = sorted(chord.project(Point(xy)) for xy in coordinates)
            if (hi-lo > tol and edge.distance(chord.interpolate(lo)) <= tol
                    and edge.distance(chord.interpolate(hi)) <= tol):
                owners.add(face)
                break
    return sorted(owners)


class _PhysicalBoundary:
    # Keep each physical ring and its precision buffer separately indexed so even components
    # meeting at a single vertex retain their boundaries.
    def __init__(self, domain):
        # Buffer individual rings: combining point-touching rings in one tiny
        # GEOS buffer can erase entire components at projected coordinates.
        rings = line_parts(domain.boundary)
        epsilon = _coordinate_tolerance(domain)
        self.index = STRtree([*rings, *(ring.buffer(epsilon) for ring in rings)])

    # Test the whole candidate edge against individual boundary rings and their buffers. Two
    # endpoints on the physical boundary can still span a cavity.
    def covers(self, edge):
        return self.index.query(edge, predicate='covered_by').size > 0


def split_edge_chain(a, b, adjacent, epsilon=None):
    """Follow a unique monotone mesh-edge chain between exact input endpoints."""
    if a not in adjacent or b not in adjacent:
        return None
    dx, dy = b[0]-a[0], b[1]-a[1]
    length = math.hypot(dx, dy)
    if epsilon is None:
        epsilon = np.spacing(max(1., *map(abs, a), *map(abs, b)))
    chain = [a]
    position = 0.
    while chain[-1] != b:
        steps = []
        for xy in adjacent[chain[-1]]:
            x, y = xy[0]-a[0], xy[1]-a[1]
            along = (x*dx + y*dy) / length
            if xy == b or (position < along < length and abs(dx*y - dy*x) <= epsilon*length):
                steps.append((xy, along))
        # A fork inside the roundoff allowance is ambiguous. Following the
        # closest branch could assign an edge to the wrong input segment.
        if len(steps) != 1:
            return None
        xy, position = steps[0]
        chain.append(xy)
    return chain


# Each clipped stream segment must survive as a continuous chain of mesh edges. Where CGAL
# supplies constraint labels, require the chain to follow those edges.
def _validate_stream_edges(vertices, triangles, neighbours, domain, constraints, constrained_edges):
    lines = [part for line in constraints for part in line_parts(line.intersection(domain))]
    if not lines:
        return
    # Limit the native path to edges explicitly marked constrained by CGAL.
    # For custom meshers only the output triangle connectivity is available.
    if constrained_edges is None:
        selected = None
    else:
        declared = np.asarray(constrained_edges, dtype=np.int64).reshape(-1, 2)
        if np.any(declared < 0) or np.any(declared >= len(vertices)):
            raise PatchError('Invalid vertex in a declared physical constraint')
        selected = {tuple(sorted(edge)) for edge in declared}
        endpoints = np.zeros(len(vertices), dtype=bool)
        endpoints[declared] = True
    adjacent, vertex_ids = {}, {}
    face_ids = np.arange(len(triangles))
    for slot in range(3):
        unique = (neighbours[:, slot] < 0) | (neighbours[:, slot] > face_ids)
        left, right = (slot+1) % 3, (slot+2) % 3
        if selected is not None:
            unique &= endpoints[triangles[:, left]] & endpoints[triangles[:, right]]
        for i, j in triangles[unique][:, [left, right]]:
            if selected is not None:
                edge = (i, j) if i < j else (j, i)
                if edge not in selected:
                    continue
                selected.remove(edge)
            a, b = tuple(vertices[i, :2]), tuple(vertices[j, :2])
            if vertex_ids.setdefault(a, i) != i or vertex_ids.setdefault(b, j) != j:
                raise PatchError('Distinct stream graph vertices have coincident coordinates')
            adjacent.setdefault(a, set()).add(b)
            adjacent.setdefault(b, set()).add(a)
    if selected:
        raise PatchError('A declared physical constraint edge is missing from the output triangles')
    epsilon = _coordinate_tolerance(domain)
    exact_edges = None
    for line in lines:
        for a, b in zip(line.coords[:-1], line.coords[1:]):
            if a == b or split_edge_chain(a[:2], b[:2], adjacent, epsilon) is not None:
                continue
            # An external mesh may represent a clipped segment as part of a
            # longer straight edge, with no vertex at the clip point.
            if constrained_edges is None and (a[:2] not in adjacent or b[:2] not in adjacent):
                if exact_edges is None:
                    exact_edges = union_all([polygon.boundary for polygon in polygons(vertices[triangles, :2])])
                if exact_edges.covers(LineString([a, b])):
                    continue
            raise PatchError(f'A physical constraint is missing from the replacement mesh: {a} to {b}')


def validate_patch(vertices, triangles, domain, physical_constraints=(), constrained_edges=None):
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
        _validate_stream_edges(vertices, triangles, neighbours, domain, physical_constraints, constrained_edges)
    return neighbours


# Translate the shared refinement settings into backend arguments so tile generation and seam
# repair use the same area, terrain and parameter criteria.
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


def _precision_interface_ids(linework, fixed_edges, epsilon):
    """Select interfaces at unresolved point/segment gaps, without snapping.

    GEOS can represent a valid cavity whose retained tip almost touches another
    constraint, but CGAL's inexact refinement constructions cannot resolve the
    slit. Only the owners of implicated interfaces should be removed. Exact
    incidences and gaps larger than coordinate roundoff are left alone.
    Check both directions: a physical vertex can approach a retained edge too.
    Use the original linework, before noding adds intersection coordinates.
    """
    if not fixed_edges:
        return np.empty(0, dtype=np.int64)
    coordinates = np.asarray([(a[:2], b[:2]) for geometry in linework for line in line_parts(geometry)
                              for a, b in zip(line.coords[:-1], line.coords[1:])])
    if not len(coordinates):
        return np.empty(0, dtype=np.int64)
    segments = linestrings(coordinates)
    fixed_coordinates, owners = [], []
    for owner, line in enumerate(fixed_edges):
        for pair in zip(line.coords[:-1], line.coords[1:]):
            fixed_coordinates.append(pair)
            owners.append(owner)
    fixed_coordinates = np.asarray(fixed_coordinates)
    owners = np.asarray(owners, dtype=np.int64)
    if not len(owners) or not len(segments):
        return np.empty(0, dtype=np.int64)

    # Use the spatial index to narrow the search, then distinguish genuine incidence from
    # distinct geometry closer than coordinate precision.
    def near_pairs(vertices, edges):
        xy = get_coordinates(vertices)
        # Bounding-box query plus a distance filter also supports Shapely 2
        # installations built with GEOS < 3.10 (which lack STRtree dwithin).
        pairs = STRtree(edges).query(boxes(xy[:, 0]-epsilon, xy[:, 1]-epsilon,
                                          xy[:, 0]+epsilon, xy[:, 1]+epsilon))
        # Distance can itself round to zero on a sloping edge. Use the exact
        # topological incidence test rather than treating distance == 0 as safe.
        selected_edges, selected_vertices = edges[pairs[1]], vertices[pairs[0]]
        return pairs[:, (~covers(selected_edges, selected_vertices)
                         & (distance(selected_edges, selected_vertices) <= epsilon))]

    fixed_vertices = points(fixed_coordinates.reshape(-1, 2))
    forward = near_pairs(fixed_vertices, segments)
    vertices = points(np.unique(coordinates.reshape(-1, 2), axis=0))
    reverse = near_pairs(vertices, linestrings(fixed_coordinates))
    return np.unique(np.concatenate((np.repeat(owners, 2)[forward[0]], owners[reverse[1]])))


def remesh_domain(domain, constraints, options, directory, fixed_edges=(), *, retained_edges=None,
                  disk_coverage=None, return_constrained_edges=False):
    """Invoke the real C++ refinement criteria, including raster constraints."""
    started = time.perf_counter()
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
    stats = {}
    try:
        (work / 'domain.json').write_text(json.dumps(mapping(domain)))
        (work / 'constraints.json').write_text(json.dumps({'type': 'FeatureCollection', 'features': features}))
        (work / 'interfaces.json').write_text(json.dumps({'edges': [list(line.coords) for line in fixed_edges],
            'export_constrained_edges': bool(constraints)}))
        output = work / 'mesh.json'
        command = _mesher_arguments(options) + ['--domain-file', str(work / 'domain.json'),
                   '--interior-plgs-file', str(work / 'constraints.json'), '--patch-output', str(output),
                   '--patch-max-vertices', str(options.get('remesh_max_vertices', 500000)),
                   '--patch-interfaces', str(work / 'interfaces.json')]
        (work / 'invocation.json').write_text(json.dumps({'command': command, 'cwd': os.getcwd(),
            'tile': options.get('tile_prefix'), 'merge': options.get('out_prefix'),
            'core_bbox': options.get('core_bbox', options.get('tile_bbox')),
            'domain_bounds': list(domain.bounds), 'domain_area': domain.area}, indent=2))
        tiny = _unbisectable_edges(noded)
        chords = _roundoff_chords(domain)
        if tiny:
            # Inspect this potentially large boundary iterator once, selecting
            # whole owners for either kind of precision defect in this cavity.
            witnesses = [*tiny, *chords]
            affected = _roundoff_owners(witnesses, retained_edges) if retained_edges is not None else []
            retained_edges = None
            on_disk = False
            if not affected and disk_coverage is not None and not disk_coverage.is_empty:
                rim = disk_coverage.boundary
                on_disk = any(rim.distance(edge) <= 2*np.spacing(max(map(abs, edge.bounds))) for edge in tiny)
            if affected or on_disk:
                error = PatchBoundaryRoundoff(
                    f'Constraint midpoint rounds to an endpoint at {list(tiny[0].coords)}; '
                    f'{len(affected)} matching local faces', witnesses, affected, preflight_only=True)
                stats['preparation_seconds'] = time.perf_counter() - started
                (work / 'backend.log').write_text(f'Backend not started: {error}\n')
                raise error
        epsilon = _coordinate_tolerance(domain)
        precision_ids = _precision_interface_ids(linework, fixed_edges, epsilon)
        if len(precision_ids):
            error = PatchPrecisionError(precision_ids, epsilon)
            stats['preparation_seconds'] = time.perf_counter() - started
            (work / 'precision.json').write_text(json.dumps({
                'interface_ids': precision_ids.tolist(), 'coordinate_tolerance': epsilon}))
            (work / 'backend.log').write_text(f'Backend not started: {error}\n')
            raise error
        if chords:
            affected = _roundoff_owners(chords, retained_edges) if retained_edges is not None else []
            if not affected:
                chords = _roundoff_chords(domain, acute_only=True)
            if chords:
                detail = (f'{len(affected)} matching retained faces' if affected else
                          'no matching local faces and a boundary angle below the mesh shape bound')
                error = PatchBoundaryRoundoff(f'{len(chords)} roundoff-width cavity components: {detail}',
                                             chords, affected, preflight_only=True)
                stats['preparation_seconds'] = time.perf_counter() - started
                (work / 'backend.log').write_text(f'Backend not started: {error}\n')
                raise error
        preparation_seconds = time.perf_counter() - started
        code, backend_seconds = _run_backend(command, work, options)
        stats_path = Path(str(output) + '.stats.json')
        stats = json.loads(stats_path.read_text()) if stats_path.exists() else {}
        stats.update(backend_seconds=backend_seconds, preparation_seconds=preparation_seconds)
        if code:
            with (work / 'backend.log').open('rb') as log:
                log.seek(0, os.SEEK_END)
                log.seek(max(0, log.tell() - 3000))
                tail = log.read().decode(errors='replace')
            message = f'C++ patch refinement failed (exit {code}): ' + tail
            quality_path = Path(str(output) + '.quality.json')
            if quality_path.exists():
                raise PatchQualityError(message, json.loads(quality_path.read_text())['faces'])
            edge_report = Path(str(output) + '.edge.json')
            if edge_report.exists():
                edge = LineString(json.loads(edge_report.read_text())['edge'])
                message = f'Constrained edge refinement failed at {list(edge.coords)}. {message}'
                if tiny or chords:
                    raise PatchBoundaryRoundoff(message, [edge, *tiny, *chords])
                raise PatchEdgeError(message, edge)
            if 'n == zone.fh' in tail and 'Refine_edges.h' in tail:
                if tiny:
                    message = f'Constraint midpoint rounds to an endpoint at {list(tiny[0].coords)}. {message}'
                chords = [*tiny, *_roundoff_chords(domain)]
                if chords:
                    raise PatchBoundaryRoundoff(message, chords)
            raise PatchError(message)
        io_started = time.perf_counter()
        data = json.loads(output.read_text())
        vertices = np.asarray(data['vertices'], dtype=float)
        triangles = np.asarray(data['triangles'], dtype=np.int64)
        io_seconds = time.perf_counter() - io_started
        validation_started = time.perf_counter()
        try:
            validate_patch(vertices, triangles, domain, constraints, data.get('constrained_edges'))
        finally:
            stats.update(io_seconds=io_seconds, validation_seconds=time.perf_counter()-validation_started)
        stats.update(data.get('stats', {}))
    except Exception as exc:
        context = options.get('tile_prefix') or options.get('out_prefix') or 'cavity'
        message = f'{context}: {exc}\nReproduction files and full backend log: {work}'
        if isinstance(exc, PatchError):
            exc.stats = dict(stats, **getattr(exc, 'stats', {}))
            exc.args = (message,)
            raise
        error = PatchError(message)
        error.stats = stats
        raise error from exc
    shutil.rmtree(work)
    return PatchResult(vertices, triangles, data.get('interfaces', []), stats,
                       data.get('constrained_edges') if return_constrained_edges else None)


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
    """Select whole input faces the tile owns and can keep outside the seam band."""
    pts = vertices[triangles, :2]
    inside = ((pts >= np.asarray(bounds[:2])) & (pts <= np.asarray(bounds[2:]))).all(axis=(1, 2))
    # A face whose centroid is outside the band can still reach across the
    # seam. Releasing the whole face leaves existing edges for the patch to meet.
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
    epsilon = _coordinate_tolerance(domain)
    affected = set()
    for failure in failures:
        if 'input_interface_ids' in failure:
            affected.update(owners[i] for i in failure['input_interface_ids'])
            continue
        triangle = Polygon(failure['vertices'])
        x0, y0, x1, y1 = triangle.bounds
        candidates = tree.query(box(x0-epsilon, y0-epsilon, x1+epsilon, y1+epsilon))
        for edge in candidates:
            if fixed_edges[edge].distance(triangle) <= epsilon:
                affected.add(owners[edge])
    return np.asarray(sorted(affected), dtype=np.int64)


def _physical_sector(failures, source_domain, constraints):
    """Certify both rays against source geometry before ruling out seam repair."""
    perimeter = None
    for face in failures:
        if 'input_rays' not in face:
            continue
        if perimeter is None:
            perimeter = source_domain.boundary
        tip = face['vertices'][face['constrained_corner']]
        rays = [LineString([tip, end]) for end in face['input_rays']]
        if all(perimeter.covers(ray) or any(line.covers(ray) for line in constraints) for ray in rays):
            return tip, face['input_sector_degrees']
    return None


def _retained_neighbours(neighbours, keep):
    """Return adjacency in the face numbering of triangles[keep].

    Filtering faces preserves their vertex order and all three edge slots.
    A surviving neighbour gets its new row number. An edge exposed by a
    removed face gets the same -1 boundary marker as an exterior edge.
    """
    if neighbours is None:
        return np.empty((0, 3), dtype=np.int64)
    mapping = np.full(len(keep), -1, dtype=np.int64)
    mapping[keep] = np.arange(np.count_nonzero(keep))
    selected = neighbours[keep]
    return np.where(selected >= 0, mapping[np.maximum(selected, 0)], -1)


def _grow_cavity(keep, neighbours, affected, rings):
    """Clear faces in keep in place and return how many were newly released.

    affected and neighbours use the face numbering of the full keep mask.
    Zero rings releases just the seeds. Following only faces still marked
    for retention prevents growth across a cavity that is already open.
    """
    frontier = np.unique(affected)
    frontier = frontier[keep[frontier]]
    removed = 0
    for _ in range(rings + 1):
        if not len(frontier):
            break
        keep[frontier] = False
        removed += len(frontier)
        adjacent = np.unique(neighbours[frontier])
        adjacent = adjacent[adjacent >= 0]
        frontier = adjacent[keep[adjacent]]
    return removed


def subdivide_parents(records, options, directory):
    """One C++ call for all split parents; failed parents receive no children."""
    work = Path(tempfile.mkdtemp(prefix='repair-', dir=directory)).resolve()
    stats = {}
    try:
        (work / 'parents.json').write_text(json.dumps({'repairs': records}))
        output = work / 'mesh.json'
        command = _mesher_arguments(options) + ['--patch-repair-file', str(work / 'parents.json'),
            '--patch-output', str(output), '--patch-max-vertices', str(options.get('remesh_max_vertices', 500000))]
        (work / 'invocation.json').write_text(json.dumps({'command': command, 'cwd': os.getcwd(),
                                                       'merge': options.get('out_prefix')}, indent=2))
        code, seconds = _run_backend(command, work, options)
        stats_path = Path(str(output) + '.stats.json')
        stats = json.loads(stats_path.read_text()) if stats_path.exists() else {}
        stats['backend_seconds'] = seconds
        if code:
            raise PatchError(f'C++ subdivision failed (exit {code}): ' + (work / 'backend.log').read_text()[-3000:])
        data = json.loads(output.read_text())
        stats.update(data.get('stats', {}))
        data['stats'] = stats
        validation_started = time.perf_counter()
        expected = {record['id']: record for record in records}
        if len(data['repairs']) != len(expected) or {r['id'] for r in data['repairs']} != set(expected):
            raise PatchError('Subdivision output lost parent identity')
        for result in data['repairs']:
            if not result['accepted']:
                continue
            record = expected[result['id']]
            points = np.asarray(result['triangles'], dtype=float).reshape((-1, 2))
            vertices, inverse = np.unique(points, axis=0, return_inverse=True)
            vertices = np.column_stack((vertices, np.zeros(len(vertices))))
            triangles = inverse.reshape((-1, 3))
            neighbours = validate_patch(vertices, triangles, Polygon(record['parent']))
            faces, slots = np.nonzero(neighbours < 0)
            actual = {tuple(sorted(map(tuple, vertices[triangles[f, [(s+1)%3, (s+2)%3]], :2])))
                      for f, s in zip(faces, slots)}
            ring = record['boundary']
            boundary = {tuple(sorted((tuple(a), tuple(b)))) for a, b in zip(ring, ring[1:]+ring[:1])}
            if actual != boundary:
                raise PatchError('Subdivision changed its prescribed boundary')
        stats['validation_seconds'] = time.perf_counter()-validation_started
    except Exception as exc:
        if isinstance(exc, PatchError):
            exc.stats = dict(stats, **getattr(exc, 'stats', {}))
            exc.args = (f'{exc}\nSubdivision reproduction files: {work}',)
            raise
        raise PatchError(f'{exc}\nSubdivision reproduction files: {work}') from exc
    shutil.rmtree(work)
    return data


def _join_exact(retained_v, retained_t, patch_v, patch_t):
    """Give the retained mesh and the replacement patch a shared vertex table.

    Each input triangle contains three row indices into its own XYZ array.
    Return (vertices, kept, added), with both sets of triangles indexing
    vertices. Each keeps its input face order so callers can still identify
    a retained parent or a generated face after the vertex IDs change.
    _stitch checks whether the remapped edges actually meet across the seam.
    """
    # A released face can leave unused vertices in retained_v. Drop those
    # before matching XY so a discarded vertex cannot supply the patch's Z.
    used = np.unique(retained_t)
    vertices = retained_v[used].tolist()
    mapping_old = np.full(len(retained_v), -1, dtype=np.int64)
    mapping_old[used] = np.arange(len(used))
    kept = mapping_old[retained_t]
    # Connectivity is defined in XY. An exact match inherits the retained
    # vertex, including its elevation. The native writer uses max_digits10
    # so shared coordinates round-trip. Snapping nearby vertices could close
    # a narrow gap, and two retained IDs at one XY make the choice ambiguous.
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
    return np.asarray(vertices, dtype=float).reshape((-1, 3)), kept, remap[patch_t]


def _repair_candidate(vertices, retained, provenance, patch, fixed, owners, options, directory, generation):
    """Stage valid children and patch interiors, removing failed two-sided regions.

    Each patch interface ID indexes fixed, and owners maps it to a row of
    retained. Its chain contains patch vertex IDs in order along that edge.
    Return None if no interface was split. Otherwise the candidate keep mask
    and neighbours index result_t. Failed parent IDs still index retained so
    the caller can release the corresponding faces in the previous mesh.
    """
    patch_v, patch_t = patch
    if len(patch.interfaces) != len(fixed):
        raise PatchError('Backend omitted retained interface lineage; rebuild meshercxx')
    patch_neighbours = build_neighbours(patch_t, len(patch_v), patch_v)
    pf, ps = np.nonzero(patch_neighbours < 0)
    boundary_owners = {tuple(sorted(patch_t[f, [(s+1)%3, (s+2)%3]])): int(f) for f, s in zip(pf, ps)}
    split, adjacent, seen = {}, {}, set()
    for interface in patch.interfaces:
        edge_id = interface['id']
        if edge_id in seen or not 0 <= edge_id < len(fixed):
            raise PatchError('Invalid interface identity')
        seen.add(edge_id)
        chain = np.asarray(interface['vertices'], dtype=np.int64)
        if len(chain) < 2 or np.any(chain < 0) or np.any(chain >= len(patch_v)):
            raise PatchError('Invalid interface subdivision')
        coords = patch_v[chain, :2]
        original = np.asarray(fixed[edge_id].coords)[[0, -1], :2]
        if not np.array_equal(coords[[0, -1]], original):
            raise PatchError('Interface endpoints changed')
        face_ids = []
        for a, b in zip(chain[:-1], chain[1:]):
            key = tuple(sorted((a, b)))
            if key not in boundary_owners:
                raise PatchError('Interface subdivision is not a patch boundary')
            face_ids.append(boundary_owners[key])
        if len(chain) > 2:
            owner = int(owners[edge_id])
            split.setdefault(owner, {})[tuple(sorted(map(tuple, original)))] = coords.tolist()
            adjacent.setdefault(owner, set()).update(face_ids)
    if not split:
        return None
    records = []
    for owner, edges in split.items():
        parent = vertices[retained[owner], :2].tolist()
        boundary = []
        for a, b in zip(parent, parent[1:]+parent[:1]):
            chain = edges.get(tuple(sorted((tuple(a), tuple(b)))), [a, b])
            if chain[0] != a:
                chain = chain[::-1]
            boundary.extend(chain[:-1])
        records.append({'id': owner, 'parent': parent, 'boundary': boundary})
    repairs = subdivide_parents(records, options, directory)
    combined_v, old_t, new_t = _join_exact(vertices, retained, patch_v, patch_t)
    ids = {tuple(v[:2]): i for i, v in enumerate(combined_v)}
    child_t, accepted, failed = [], [], []
    for repair in repairs['repairs']:
        owner = repair['id']
        if repair['accepted']:
            accepted.append(owner)
            for triangle in repair['triangles']:
                # Children may only use supplied boundary vertices.
                child_t.append([ids[tuple(point)] for point in triangle])
        else:
            failed.append(owner)
    selection = np.ones(len(old_t) + len(new_t), dtype=bool)
    selection[accepted] = False
    all_t = np.vstack((old_t, new_t))
    # Removing accepted parents shifts the rows of surviving faces. Translate
    # the failed parents and their patch-side neighbours into result_t before
    # opening a cavity on both sides of each failed subdivision.
    old_to_new = np.cumsum(selection) - 1
    result_t = np.vstack((all_t[selection], np.asarray(child_t, dtype=np.int64).reshape((-1, 3))))
    sources = np.concatenate((provenance, np.full(len(new_t), generation)))[selection]
    sources = np.concatenate((sources, np.full(len(child_t), generation)))
    keep = np.ones(len(result_t), dtype=bool)
    neighbours = build_neighbours(result_t, len(combined_v), combined_v)
    if failed:
        seeds = list(old_to_new[failed])
        for owner in failed:
            seeds.extend(old_to_new[len(old_t) + np.asarray(sorted(adjacent[owner]), dtype=np.int64)])
        _grow_cavity(keep, neighbours, seeds, 1)
    return combined_v, result_t, sources, keep, neighbours, repairs['stats'], len(accepted), failed


def _retained_boundary_edges(vertices, triangles, adjacency):
    """Yield each exposed edge's owner row in triangles and its two XY endpoints.

    Neighbour slot i lies opposite vertex i. When triangles is a compact
    subset, the caller maps the owner back through retained_ids before
    releasing it in the full working mesh.
    """
    for face, slot in zip(*np.nonzero(adjacency < 0)):
        yield face, vertices[triangles[face, [(slot+1)%3, (slot+2)%3]], :2]


def _conforming_remainder(vertices, triangles, neighbours, keep):
    """Shrink keep in place until the surviving boundary conforms, then return its footprint.

    A release exposes new edges. Repeating the audit catches hanging nodes
    that were interior to the previous selection before it becomes the input
    to the next repair.
    """
    while True:
        face_ids = np.flatnonzero(keep)
        remaining = triangles[face_ids]
        adjacency = _retained_neighbours(neighbours, keep)
        coincident, hanging = boundary_node_conflicts(vertices, remaining, adjacency, 1e-8)
        nodes = set(coincident)
        coarse_faces = []
        for face, _, interior_nodes in hanging:
            coarse_faces.append(face)
            nodes.update(interior_nodes)
        if not nodes:
            return mesh_coverage(vertices, remaining, adjacency)
        # IDs from the audit index the compact remainder. Remove the edge's
        # owner and every surviving triangle incident to an offending node.
        remove = np.isin(remaining, list(nodes)).any(axis=1)
        remove[coarse_faces] = True
        keep[face_ids[remove]] = False


def _stitch(retained_v, retained_t, patch_v, patch_t, domain, physical_boundary=None,
            retained_neighbours=None, result_neighbours=None):
    """Join the vertex tables and check that the faces meet across the cavity.

    Kept faces occupy the first rows of the result. An exposed retained edge
    inside the domain needs a neighbour in the added rows. A split edge can
    share both endpoints yet still leave a hanging node between them.
    InterfaceChanged reports rows of retained_t for the caller to reopen.
    """
    result_v, kept, added = _join_exact(retained_v, retained_t, patch_v, patch_t)
    result_t = np.vstack((kept, added))
    neighbours = (build_neighbours(result_t, len(result_v), result_v)
                  if result_neighbours is None else result_neighbours)
    if len(kept):
        if physical_boundary is None:
            physical_boundary = _PhysicalBoundary(domain)
        old_neighbours = (build_neighbours(kept, len(result_v), result_v)
                          if retained_neighbours is None else retained_neighbours)
        faces, slots = np.nonzero(old_neighbours < 0)
        changed = []
        for face, slot in zip(faces, slots):
            edge = result_v[kept[face, [(slot + 1) % 3, (slot + 2) % 3]], :2]
            line = LineString(edge)
            if not physical_boundary.covers(line) and neighbours[face, slot] < len(kept):
                changed.append(face)
        if changed:
            raise InterfaceChanged(changed)
    coincident, hanging = boundary_subdivisions(result_v, result_t, neighbours, 1e-8)
    if coincident or next(hanging, None) is not None:
        raise PatchError('Joined mesh has coincident or hanging boundary nodes')
    return result_v, result_t


def merge_cavity(vertices_a, triangles_a, bounds_a, vertices_b, triangles_b, bounds_b,
                 domain, constraints, options, directory, mesher=remesh_domain,
                 *, prepared_domain=None, initial_keep=None, disk_coverage=None, phase_prefix=''):
    """Repair a private working mesh; publish only a conforming complete result."""
    started = time.perf_counter()
    axis, seam = _shared_interface(bounds_a, bounds_b)
    bounds = [min(bounds_a[0], bounds_b[0]), min(bounds_a[1], bounds_b[1]),
              max(bounds_a[2], bounds_b[2]), max(bounds_a[3], bounds_b[3])]
    supplied_width = options.get('remesh_band_width')
    width = float(2 * math.sqrt(options['max_area']) if supplied_width is None else supplied_width)
    attempts = int(options.get('remesh_max_attempts', 12))
    if not np.isfinite(width) or width <= 0 or attempts < 1:
        raise PatchError('Remesh band width and attempt count must be positive')
    context = options.get('out_prefix') or options.get('tile_prefix')
    source_domain = domain
    core_area = domain.intersection(box(*bounds)).area
    # Intermediate group boxes can recreate the same acute cuts as initial
    # tiles. Give only those corners local halos; ownership boxes never change.
    domain = (prepared_domain if prepared_domain is not None else
              _clipped_domain(domain, bounds, 0, width/2, context, options.get('mesher_debug', False)))
    if domain.is_empty:
        return np.empty((0, 3)), np.empty((0, 3), dtype=np.int64), {'attempts': 0}
    options = dict(options, core_bbox=bounds)
    joined_v = np.vstack((vertices_a, vertices_b))
    joined_t = np.vstack((triangles_a, triangles_b + len(vertices_a)))
    sources = np.arange(len(joined_t))  # Nonnegative: unchanged original; negative: generation.
    keep = (np.asarray(initial_keep, dtype=bool).copy() if initial_keep is not None else
            np.concatenate((_retain_mask(vertices_a, triangles_a, bounds_a, axis, seam, width),
                            _retain_mask(vertices_b, triangles_b, bounds_b, axis, seam, width))))
    if keep.shape != (len(joined_t),):
        raise PatchError('Initial cavity selection does not match the working faces')
    neighbours = build_neighbours(joined_t, len(joined_v), joined_v) if len(joined_t) else None
    physical_boundary = _PhysicalBoundary(domain)
    growth_rings, interface_recovery = 1, False
    label = f' [{Path(context).name}]' if context else ''
    last_error = None
    history = []
    failed_parent_shapes = set()
    grow_failed_parents = False
    lineage_available = mesher is remesh_domain
    totals = {'subdivided_parents': 0, 'local_repairs': 0, 'patch_triangles_generated': 0,
              'failed_subdivision_parents': 0, 'candidate_interface_faces': 0}

    # Release connected retained components along artificial interfaces on the last recovery.
    # Components bounded only by physical edges can stay intact.
    def final_recovery():
        nonlocal interface_recovery
        retained_ids = np.flatnonzero(keep)
        retained = joined_t[keep]
        if not len(retained):
            return
        adjacency = _retained_neighbours(neighbours, keep)
        faces, slots = np.nonzero(adjacency < 0)
        owners = [f for f, s in zip(faces, slots)
                  if not physical_boundary.covers(LineString(joined_v[retained[f, [(s+1)%3, (s+2)%3]], :2]))]
        if owners:
            removed = _grow_cavity(keep, neighbours, retained_ids[np.unique(owners)], len(joined_t))
            interface_recovery = True
            debug_log(options, f'Final interface recovery{label}: removed={removed}, retained={np.count_nonzero(keep)}; '
                      'remeshing components touching the cavity with unchanged quality checks and vertex budget')

    # Use face provenance to distinguish original triangles from reused patch work and
    # accumulate costs from every attempt into the merge report.
    def finish(v, t, provenance, attempt, patch_count):
        stats = dict(totals, attempts=attempt+1, band_width=width, patch_triangles=patch_count,
                     boundary_halo_area=domain.area-core_area,
                     retained_triangles=int(np.count_nonzero(provenance > -(attempt+1))),
                     original_triangles_retained=int(np.count_nonzero(provenance >= 0)),
                     generated_triangles_reused=int(np.count_nonzero((provenance < 0) & (provenance > -(attempt+1)))),
                     interface_recovery=interface_recovery, attempt_stats=history,
                     elapsed_seconds=time.perf_counter()-started)
        for key in ('backend_seconds', 'preparation_seconds', 'validation_seconds', 'stitch_seconds',
                    'io_seconds', 'repair_seconds', 'refinement_seconds', 'refinement_steps',
                    'quality_repair_vertices', 'raster_seconds', 'raster_evaluations'):
            stats[key] = sum(record.get(key, 0) for record in history)
        debug_log(options, f'Merge complete{label}: original_retained={stats["original_triangles_retained"]}, '
                  f'generated_reused={stats["generated_triangles_reused"]}, '
                  f'subdivided={stats["subdivided_parents"]}, elapsed={stats["elapsed_seconds"]:.2f}s')
        return v, t, stats

    attempts_used = 0
    while attempts_used < attempts:
        attempt = attempts_used
        if prepared_domain is not None and len(joined_t) > options.get('remesh_max_active_faces', 2000000):
            raise PatchError('Seam active-face budget exceeded during local repair; '
                             'increase mpi_remesh_max_active_faces only if memory permits')
        if len(history) == attempt:
            if attempt > 0 and attempt+1 == attempts and not interface_recovery:
                final_recovery()
            history.append({'attempt': attempt+1})
        record = history[attempt]
        prepare_started = time.perf_counter()
        retained_ids = np.flatnonzero(keep)
        # An outer disk-interface attempt may need several local cavity solves
        # before it can offer the disk mesh a matching edge.
        report_phase(options, f'{phase_prefix}solving cavity, attempt {attempt+1}/{attempts} '
                     f'({len(retained_ids):,} retained of {len(joined_t):,} local faces)')
        retained, retained_sources = joined_t[keep], sources[keep]
        adjacency = _retained_neighbours(neighbours, keep)
        coverage = mesh_coverage(joined_v, retained, adjacency)
        if coverage.difference(domain).area > area_tolerance(domain):
            raise PatchError('Retained triangles extend outside the intended domain')
        cavity = polygon_domain(domain.difference(coverage))
        record.update(cavity_area=cavity.area, retained_triangles=len(retained))
        debug_log(options, f'Remesh attempt {attempt+1}{label}: width={width:g}, retained={len(retained)}, cavity_area={cavity.area:g}')
        fixed, owners = [], []
        cavity_boundary = cavity.boundary
        for f, coords in _retained_boundary_edges(joined_v, retained, adjacency):
            line = LineString(coords)
            if not physical_boundary.covers(line) and cavity_boundary.covers(line):
                fixed.append(line)
                owners.append(f)
        record['preparation_seconds'] = record.get('preparation_seconds', 0) + (time.perf_counter()-prepare_started)
        attempts_used += 1
        try:
            if mesher is remesh_domain:
                boundary = _retained_boundary_edges(joined_v, retained, adjacency)
                patch = None
                # With zero retained faces there is no retained edge lineage
                # to translate. Partition only this expensive rebuild, leaving
                # the ordinary local repair and final disk checks in control.
                if (not len(retained_ids) and '_cavity_team' in options
                        and len(joined_t) > options['cavity_target_faces'] and not cavity.is_empty):
                    from .parallel_cavity import solve
                    patch = solve(cavity, constraints, options, 1-axis, len(joined_t), phase_prefix.rstrip(' |'))
                    if patch is None:
                        report_phase(options, f'{phase_prefix}solving unsplit cavity ({len(joined_t):,} local faces)')
                if patch is None:
                    patch = mesher(cavity, constraints, options, directory, fixed,
                                   retained_edges=boundary, disk_coverage=disk_coverage)
            else:
                patch = mesher(cavity, constraints, options, directory, fixed)
            patch_v, patch_t = patch
            if isinstance(patch, PatchResult):
                lineage_available = True
                for key, value in patch.stats.items():
                    record[key] = record.get(key, 0) + value
            elif len(patch_t):
                validation_started = time.perf_counter()
                validate_patch(patch_v, patch_t, cavity, constraints)
                record['validation_seconds'] = time.perf_counter()-validation_started
            totals['patch_triangles_generated'] += len(patch_t)
            repair = None
            if isinstance(patch, PatchResult) and fixed:
                repair_started = time.perf_counter()
                try:
                    repair = _repair_candidate(joined_v, retained, retained_sources, patch, fixed, owners,
                                               options, directory, -(attempt+1))
                finally:
                    record['repair_seconds'] = time.perf_counter()-repair_started
            stitch_started = time.perf_counter()
            if repair is not None:
                v, t, provenance, selection, candidate_neighbours, repair_stats, accepted, failed_owners = repair
                failed = len(failed_owners)
                for key, value in repair_stats.items():
                    record[key] = record.get(key, 0) + value
                changed_faces = np.empty(0, dtype=np.int64)
                if not failed:
                    try:
                        output_v, output_t = _stitch(v, t, np.empty((0, 3)), np.empty((0, 3), dtype=np.int64),
                            domain, physical_boundary, candidate_neighbours, candidate_neighbours)
                    except InterfaceChanged as exc:
                        # These IDs belong to t, which includes generated faces
                        # and replaces subdivided parents. They MUST NOT reach
                        # the outer handler, whose retained_ids index the earlier
                        # working mesh. Remove only affected candidate faces and
                        # their neighbours, then validate before promotion.
                        changed_faces = exc.faces
                        if attempt+1 == attempts:
                            raise PatchError(f'{len(changed_faces)} candidate interface faces require another local repair') from exc
                        _grow_cavity(selection, candidate_neighbours, changed_faces, 1)
                    else:
                        totals['subdivided_parents'] += accepted
                        record['stitch_seconds'] = time.perf_counter()-stitch_started
                        return finish(output_v, output_t, provenance, attempt, len(patch_t))
                if attempt+1 == attempts:
                    raise PatchError(f'{failed} parent subdivisions require another local repair')
                if failed:
                    # Generated face IDs change with every candidate. An exact
                    # XY triangle identifies the same failed parent across
                    # generations even when its index or orientation changes.
                    if not grow_failed_parents:
                        shapes = {tuple(sorted(map(tuple, joined_v[retained[i], :2]))) for i in failed_owners}
                        grow_failed_parents = not failed_parent_shapes.isdisjoint(shapes)
                        failed_parent_shapes.update(shapes)
                    if grow_failed_parents:
                        # Keep the previously validated remainder, and enlarge
                        # its cavity. A failed candidate can no longer put back
                        # the parents we just released and repeat the cycle.
                        removed = _grow_cavity(keep, neighbours, retained_ids[failed_owners], growth_rings)
                        record['repeated_parent_growth_faces'] = removed
                        totals['failed_subdivision_parents'] += failed
                        growth_rings = min(2*growth_rings, 4)
                        debug_log(options, f'Growing cavity after repeated parent failure{label}: '
                                  f'failed_parents={failed}, removed={removed}')
                        continue
                # Validate the staged remainder BEFORE promoting it. Every
                # preserved generated triangle passed backend quality; all
                # children passed both quality and exact footprint checks.
                before = np.count_nonzero(selection)
                remainder = _conforming_remainder(v, t, candidate_neighbours, selection)
                boundary_faces = int(before - np.count_nonzero(selection))
                if boundary_faces:
                    record['staged_boundary_faces_released'] = boundary_faces
                    debug_log(options, f'Releasing {boundary_faces} staged faces at boundary conflicts{label}')
                if remainder.difference(domain).area > area_tolerance(domain):
                    raise PatchError('Staged repair extends beyond the domain')
                joined_v, joined_t, sources, keep = v, t, provenance, selection
                neighbours = candidate_neighbours
                totals['local_repairs'] += 1
                preserved = int(np.count_nonzero((sources < 0) & keep))
                totals['subdivided_parents'] += accepted
                totals['failed_subdivision_parents'] += failed
                totals['candidate_interface_faces'] += len(changed_faces)
                debug_log(options, f'Preserving patch work{label}: generated={preserved}, '
                          f'failed_parents={failed}, candidate_interface_faces={len(changed_faces)}, '
                          f'removed={np.count_nonzero(~keep)}')
                record.update(stitch_seconds=time.perf_counter()-stitch_started,
                              generated_preserved=preserved, failed_subdivisions=failed,
                              candidate_interface_faces=len(changed_faces))
                continue
            output_v, output_t = _stitch(joined_v, retained, patch_v, patch_t, domain, physical_boundary,
                                         retained_neighbours=adjacency)
            record['stitch_seconds'] = time.perf_counter()-stitch_started
            return finish(output_v, output_t, np.concatenate((retained_sources, np.full(len(patch_t), -(attempt+1)))),
                          attempt, len(patch_t))
        except (PatchError, MeshTopologyError) as exc:
            last_error = exc
            for key, value in getattr(exc, 'stats', {}).items():
                record[key] = record.get(key, 0) + value
            if isinstance(exc, PatchBoundaryRoundoff):
                # A subdivided physical perimeter can form a one-ULP ring
                # even when there are no fixed interfaces. Release whole rim
                # faces before spending a refinement attempt. Every pass
                # strictly shrinks keep, including on the final allowed solve.
                affected = exc.owners
                if affected is None and len(retained):
                    affected = _roundoff_owners(exc.chords, _retained_boundary_edges(joined_v, retained, adjacency))
                if affected and (exc.preflight_only or attempt+1 < attempts):
                    keep[retained_ids[affected]] = False
                    record['roundoff_parents_released'] = record.get('roundoff_parents_released', 0) + len(affected)
                    debug_log(options, f'Releasing {len(affected)} retained faces along a roundoff-width boundary{label}')
                    if exc.preflight_only:
                        attempts_used = attempt
                    continue
                exc.preflight_only = attempt == 0 and exc.preflight_only
                raise
            record['rejection'] = str(exc)
            debug_log(options, f'Remesh candidate rejected{label}: {exc}')
            if isinstance(exc, PatchQualityError):
                physical = _physical_sector(exc.failures, source_domain, constraints)
                if physical is not None:
                    tip, degrees = physical
                    raise PhysicalConstraintAngle(
                        f'Physical constraints at {tip} enclose {degrees:.6g} degrees, below the mesh shape bound. '
                        f'{exc}', exc.failures) from exc
            if (attempt+1 == attempts or not len(retained) or isinstance(exc, BackendTerminated)
                    or 'Patch vertex budget exceeded' in str(exc)):
                break
            if attempt+2 == attempts:
                continue  # The final attempt handles all remaining artificial interfaces.
            local_faces = np.empty(0, dtype=np.int64)
            if isinstance(exc, InterfaceChanged):
                local_faces = exc.faces
            elif isinstance(exc, PatchQualityError):
                local_faces = _quality_interface_faces(fixed, owners, exc.failures, domain)
            elif isinstance(exc, PatchEdgeError):
                local_faces = _roundoff_owners([exc.edge], _retained_boundary_edges(joined_v, retained, adjacency))
            elif isinstance(exc, PatchPrecisionError):
                local_faces = np.unique(np.asarray(owners)[exc.interface_ids])
            if len(local_faces):
                affected = retained_ids[local_faces]
                removed = _grow_cavity(keep, neighbours, affected, growth_rings)
                debug_log(options, f'Growing cavity locally{label}: affected={len(affected)}, '
                          f'neighbour_rings={growth_rings}, removed={removed}')
                # Once work can be preserved, avoid the previous 64/128-ring
                # erosion. Legacy/custom meshers without lineage still use the
                # original recovery path for compatibility.
                growth_rings = min(2*growth_rings, 4 if lineage_available else len(joined_t))
                continue
            growth_rings = 1
            width *= 2
            keep &= (np.abs(joined_v[joined_t, axis] - seam) >= width).all(axis=1)
    raise PatchError(f'No valid cavity replacement after {attempt+1} attempts{label}: {last_error}') from last_error
