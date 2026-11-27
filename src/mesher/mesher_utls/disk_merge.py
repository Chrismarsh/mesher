"""Merge an expandable seam region; immutable interiors remain file-backed.

Full-mesh passes are sequential/batched IO (selection, topology, publication).
Only the seam region and its interface/physical-boundary geometry are assembled
in memory. The existing native refinement and parent-repair rules still apply.
"""
from contextlib import ExitStack
import math
import json
import zlib
import os
from pathlib import Path
import tempfile
import time
import zipfile

import numpy as np
from shapely import STRtree, polygons, union_all
from shapely.geometry import LineString, Polygon, box
from shapely.ops import polygonize

from mesher.mesh_store import batches, disk_array, open_mesh_arrays, check_boundary_geometry
from mesher.mesh_topology import MeshTopologyError, build_neighbours
from . import mpi_remesh as remesh
from .mpi_progress import debug_log, report_phase
from .mpi_resume import write_marker


def selected_boundary(vertices, triangles, neighbours, active):
    """Walk the perimeter of the faces that will remain on disk.

    Yield the disk face ID, its edge's vertex IDs, their XY coordinates and
    a flag marking an active neighbour. Only flagged edges need partners
    in the local mesh. Both kinds of edge bound the disk footprint.
    """
    for lo, hi in batches(len(triangles)):
        block = neighbours[lo:hi]
        retained = ~active[lo:hi]
        for slot in range(3):
            adjacent = block[:, slot]
            cut = adjacent < 0
            valid = ~cut
            cut[valid] = active[adjacent[valid]]
            for local in np.flatnonzero(retained & cut):
                face = lo + int(local)
                tri = triangles[face]
                p, q, r = vertices[tri, :2]
                i, j = tri[(slot+1) % 3], tri[(slot+2) % 3]
                a, b = vertices[i, :2], vertices[j, :2]
                if (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0]) < 0:
                    i, j, a, b = j, i, b, a
                yield face, (i, j), (tuple(a), tuple(b)), adjacent[local] >= 0


def inactive_coverage(vertices, triangles, neighbours, active):
    """Return the disk footprint and the two lookups needed to reattach it.

    Sorted XY endpoint pairs identify the disk face behind each cut edge
    even after the local mesh is renumbered. The vertex-ID-to-XY map lets
    a later expansion recover shared vertices from this input tile.
    """
    directed, edges, interfaces = set(), [], {}
    cut_vertices = {}
    for face, ids, pair, internal in selected_boundary(vertices, triangles, neighbours, active):
        directed.add(pair)
        edges.append(LineString(pair))
        if internal:
            interfaces[tuple(sorted(pair))] = face
            cut_vertices.update(zip(ids, pair))
    pieces = []
    # polygonize also returns polygons inside holes. selected_boundary orients
    # each edge with its disk face on the left, so a CCW exterior that follows
    # one of those edges encloses retained material.
    for polygon in polygonize(edges):
        coords = list(polygon.exterior.coords)
        if not polygon.exterior.is_ccw:
            coords.reverse()
        if (coords[0], coords[1]) in directed:
            pieces.append(polygon)
    coverage = remesh.polygon_domain(union_all(pieces))
    total = 0.
    for lo, hi in batches(len(triangles)):
        pts = vertices[triangles[lo:hi][~active[lo:hi]], :2]
        u, v = pts[:, 1]-pts[:, 0], pts[:, 2]-pts[:, 0]
        total += .5*np.abs(u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]).sum()
    if abs(total-coverage.area) > remesh.area_tolerance(coverage):
        raise remesh.PatchError('Disk interior has inconsistent boundary coverage')
    return coverage, interfaces, cut_vertices


def grow_active(active, neighbours, seeds, rings, max_faces=None):
    """Expand a local frontier without allocating a whole-domain index vector."""
    frontier = np.asarray(seeds, dtype=np.int64)
    added = []
    count = 0
    for _ in range(rings+1):
        frontier = frontier[~active[frontier]]
        if not len(frontier):
            break
        count += len(frontier)
        if max_faces is not None and count > max_faces:
            # The caller is testing whether a whole component fits. Leaving a
            # partial component active would exclude those faces from the disk mesh
            # even though no face IDs are returned for the local solve.
            for ring in added:
                active[ring] = False
            return None
        active[frontier] = True
        added.append(frontier)
        adjacent = neighbours[frontier].ravel()
        frontier = np.unique(adjacent[adjacent >= 0])
    return np.concatenate(added) if added else np.empty(0, dtype=np.int64)


def close_disk_components(inputs, active, pending, capacity):
    """Finish growing affected collars that fit, leaving oversized ones intact."""
    closed = False
    for side, ((_, _, neighbours), selection, (added, removed)) in enumerate(zip(inputs, active, pending)):
        if not len(added) or capacity <= 0:
            continue
        adjacent = neighbours[added].ravel()
        seeds = np.unique(adjacent[adjacent >= 0])
        extra = grow_active(selection, neighbours, seeds, len(selection), max_faces=capacity)
        if extra is not None and len(extra):
            pending[side] = (np.concatenate((added, extra)), removed)
            capacity -= len(extra)
            closed = True
    return closed


def widen_disk_collars(inputs, active, pending, capacity, rings):
    """Add complete retained rings on both tiles while the cached solve fits."""
    frontiers = [added for added, _ in pending]
    collars = [[added] for added in frontiers]
    for _ in range(rings):
        if capacity <= 0:
            break
        next_ring = []
        for (_, _, neighbours), selection, frontier in zip(inputs, active, frontiers):
            adjacent = neighbours[frontier].ravel()
            candidates = np.unique(adjacent[adjacent >= 0])
            next_ring.append(candidates[~selection[candidates]])
        count = sum(map(len, next_ring))
        # Both tiles draw from one worker's budget. Check the combined ring
        # before changing either mask so tile order cannot favour one side.
        if count == 0 or count > capacity:
            break
        for selection, collar, faces in zip(active, collars, next_ring):
            selection[faces] = True
            collar.append(faces)
        capacity -= count
        frontiers = next_ring
    for side, collar in enumerate(collars):
        if len(collar) > 1:
            pending[side] = (np.concatenate(collar), pending[side][1])


def compact_region(vertices, triangles, active):
    """Return local vertices, locally indexed triangles and their input face IDs.

    The third array maps each local face back to its row in the tile archive.
    It lets a retry release a specific disk parent after compaction has
    renumbered the working mesh.
    """
    # This allocation is proportional to the active seam, never the interiors.
    ids = np.concatenate([np.flatnonzero(active[lo:hi])+lo for lo, hi in batches(len(triangles))]) if len(triangles) else np.empty(0, dtype=np.int64)
    selected = triangles[ids]
    used, inverse = np.unique(selected, return_inverse=True)
    return vertices[used], inverse.reshape(-1, 3), ids


def boundary_index(vertices, triangles, neighbours=None):
    """Index exposed edges by XY so independently numbered meshes can be compared.

    Return unordered endpoint pairs for matching whole edges and a map from
    endpoint XY to local vertex ID for attaching the next ring of disk faces.
    """
    if not len(triangles):
        return set(), {}
    if neighbours is None:
        neighbours = build_neighbours(triangles, len(vertices), vertices)
    faces, slots = np.nonzero(neighbours < 0)
    pairs = triangles[faces[:, None], (slots[:, None] + [1, 2]) % 3]
    points = {i: tuple(vertices[i, :2]) for i in np.unique(pairs)}
    edges = {tuple(sorted((points[a], points[b]))) for a, b in pairs}
    return edges, {xy: i for i, xy in points.items()}


def extend_cached_region(cached, inputs, pending):
    """Extend a solved seam with newly active disk faces for the next retry.

    Each tile has its own vertex IDs. Its cut_ids map lets a newly loaded
    triangle reuse the cached vertex row at each known shared endpoint.
    The returned keep mask retains the cached faces and excludes the disk
    faces marked for replacement from the newly loaded collar.
    """
    # Faces surrounding a rejected split must be free to change on both
    # sides. Keeping the seam fan would pin the same split on the next try.
    local_v, local_t, shared, release = cached
    vertices, triangles = [local_v], [local_t]
    keep = [np.ones(len(local_t), dtype=bool)]
    keep[0][release] = False
    offset = len(local_v)
    for (v, t, _), (added, removed), cut_ids in zip(inputs, pending, shared):
        used, inverse = np.unique(t[added], return_inverse=True)
        ids = np.fromiter((cut_ids.get(i, -1) for i in used), dtype=np.int64, count=len(used))
        new = ids < 0
        new_vertices = v[used[new]]
        vertices.append(new_vertices)
        ids[new] = offset + np.arange(len(new_vertices))
        offset += len(new_vertices)
        triangles.append(ids[inverse].reshape(-1, 3))
        keep.append(~np.isin(added, removed))
    return np.vstack(vertices), np.vstack(triangles), np.concatenate(keep)


def subdivide_disk_interface(inputs, active, parents, missing, vertices, triangles, boundary_ids, keys, options, directory):
    """Stage each repairable parent. Only its own complete chain can replace it."""
    adjacent = {}
    for a, b in keys:
        adjacent.setdefault(a, []).append(b)
        adjacent.setdefault(b, []).append(a)
    records, owners = [], []
    expected_faces = len(triangles)
    boundary_vertices = 0
    for side, ((v, t, _), faces, edges) in enumerate(zip(inputs, parents, missing)):
        for face in faces:
            corners = list(map(tuple, v[t[face], :2]))
            ring = []
            for a, b in zip(corners, corners[1:] + corners[:1]):
                edge = tuple(sorted((a, b)))
                if edge not in edges:
                    ring.append(a)
                    continue
                chain = remesh.split_edge_chain(a, b, adjacent)
                if chain is None:
                    break
                ring.extend(chain[:-1])
            else:
                # Each ring follows one triangle's perimeter. With no interior
                # vertices, triangulating n boundary vertices yields n-2 faces.
                if (expected_faces + len(ring)-2 > options.get('remesh_max_active_faces', 2000000)
                        or boundary_vertices + len(ring) > options.get('remesh_max_vertices', 500000)):
                    continue
                expected_faces += len(ring)-2
                boundary_vertices += len(ring)
                records.append({'id': len(records), 'parent': corners, 'boundary': ring})
                owners.append((side, face))
    if not records:
        return None

    result = remesh.subdivide_parents(records, options, directory)
    # One rejected child keeps its whole parent on disk. Other parents in
    # the same backend batch can still contribute all of their children.
    accepted = [[], []]
    ids = dict(boundary_ids)
    extra, children = [], []
    for repair in result['repairs']:
        if not repair['accepted']:
            continue
        side, face = owners[repair['id']]
        accepted[side].append(face)
        v, t, _ = inputs[side]
        for point in v[t[face]]:
            xy = tuple(point[:2])
            if xy not in ids:
                ids[xy] = len(vertices) + len(extra)
                extra.append(point)
        children.extend([ids[tuple(p)] for p in child] for child in repair['triangles'])
    if not children:
        return None
    if extra:
        vertices = np.vstack((vertices, extra))
    triangles = np.vstack((triangles, children))
    try:
        repaired_keys, repaired_ids = boundary_index(vertices, triangles)
    except MeshTopologyError:
        return None

    # Moving an accepted parent into the seam exposes its other edges and may
    # introduce a shared corner. Re-index the new cut before caching any work.
    unresolved, shared = [], []
    for values, selection, faces, old_edges in zip(inputs, active, accepted, missing):
        edges, endpoints = set(), {}
        selection[faces] = True
        try:
            for _, indices, pair, internal in selected_boundary(*values, selection):
                if internal:
                    edge = tuple(sorted(pair))
                    if edge not in repaired_keys:
                        if edge not in old_edges:
                            return None
                        edges.add(edge)
                    endpoints.update(zip(indices, pair))
        finally:
            selection[faces] = False
        unresolved.append(edges)
        shared.append(endpoints)
    return vertices, triangles, accepted, unresolved, shared, repaired_ids, repaired_keys


def _split_edge_fans(triangles, keys, boundary_ids, missing):
    """Find all faces incident to inserted vertices on a split disk edge.

    Several edges can share a triangle. Their fans must be retriangulated
    together, otherwise two accepted edits would replace the same face.
    """
    adjacent = {}
    for a, b in keys:
        adjacent.setdefault(a, []).append(b)
        adjacent.setdefault(b, []).append(a)
    chains = {}
    at_vertex = {}
    for edge in sorted(set.union(*missing)):
        chain = remesh.split_edge_chain(*edge, adjacent)
        if chain is None or len(chain) < 3:
            continue
        chains[edge] = chain
        for xy in chain[1:-1]:
            at_vertex.setdefault(boundary_ids[xy], []).append(edge)
    if not chains:
        return

    # Scan the local mesh once for all fans. The input archive interiors
    # are never indexed here, even when just one edge remains on a huge tile.
    faces_at_edge = {e: [] for e in chains}
    edges_at_face = {}
    for face in np.flatnonzero(np.isin(triangles, list(at_vertex)).any(axis=1)):
        edges = {e for v in triangles[face] for e in at_vertex.get(v, ())}
        edges_at_face[face] = edges
        for edge in edges:
            faces_at_edge[edge].append(face)

    unseen = set(chains)
    for edge in chains:
        if edge not in unseen:
            continue
        stack, group, faces = [edge], [], set()
        unseen.remove(edge)
        while stack:
            current = stack.pop()
            group.append(chains[current])
            for face in faces_at_edge[current]:
                if face not in faces:
                    faces.add(face)
                    neighbours = unseen.intersection(edges_at_face[face])
                    unseen.difference_update(neighbours)
                    stack.extend(neighbours)
        yield sorted(faces), group


def repair_seam_boundary(vertices, triangles, keys, boundary_ids, missing, cut_vertices,
                         constraints, options, directory):
    """Try the seam side of a rejected disk subdivision with its other edges fixed.

    Replacing a fan can remove the extra seam vertex without changing the
    adjacent disk triangle. Only accept native, quality-tested triangles with
    the exact required perimeter. Return face IDs still needing two-sided
    remeshing in the numbering of the returned triangle array.
    """
    required = {xy for side in cut_vertices for xy in side.values()}
    streams = STRtree(constraints) if constraints else None
    records, fans, owners = [], [], []
    retry = np.zeros(len(triangles), dtype=bool)
    budget = options.get('remesh_max_vertices', 500000)

    def perimeter(xys):
        return {tuple(sorted((a, b))) for a, b in zip(xys, xys[1:] + xys[:1])}

    for faces, chains in _split_edge_fans(triangles, keys, boundary_ids, missing):
        retry[faces] = True
        interior = {xy for chain in chains for xy in chain[1:-1]}
        if not required.isdisjoint(interior):
            continue
        original = union_all(polygons(vertices[triangles[faces], :2]))
        if original.geom_type != 'Polygon' or original.interiors:
            continue
        # Removing a vertex is legal only along the verified chains. A
        # branching boundary, hole, or physical stream needs the usual
        # constrained two-sided solve. Never snap endpoints to make them fit.
        outline = list(original.exterior.coords[:-1])
        ring = [xy for xy in outline if xy not in interior]
        if len(ring) < 3 or len(ring) > budget:
            continue
        old = perimeter(outline)
        child_edges = {tuple(sorted((a, b))) for chain in chains for a, b in zip(chain, chain[1:])}
        restored = {tuple(sorted((ch[0], ch[-1]))) for ch in chains}
        if not child_edges <= old or perimeter(ring) != ((old - child_edges) | restored):
            continue
        target = Polygon(ring)
        if not target.is_valid or original.symmetric_difference(target).area > remesh.area_tolerance(original):
            continue
        if streams is not None and (streams.query(original, predicate='intersects').size or
                                    streams.query(target, predicate='intersects').size):
            continue
        record = {'id': len(records), 'parent': outline, 'boundary': ring}
        records.append(record)
        budget -= len(ring)
        fans.append(faces)
        owners.append({tuple(vertices[v, :2]): v for v in np.unique(triangles[faces])})

    if not records:
        return vertices, triangles, missing, boundary_ids, np.flatnonzero(retry)
    children = []
    drop = np.zeros(len(triangles), dtype=bool)
    for result in remesh.subdivide_parents(records, options, directory)['repairs']:
        if result['accepted']:
            i = result['id']
            drop[fans[i]] = True
            children.extend([[owners[i][tuple(xy)] for xy in tri] for tri in result['triangles']])
    if not children:
        return vertices, triangles, missing, boundary_ids, np.flatnonzero(retry)

    # Compaction drops the removed split vertices as well. They must not
    # reappear as unused nodes in the final CHM archive or a cached retry.
    used, inverse = np.unique(np.vstack((triangles[~drop], children)), return_inverse=True)
    updated_v, updated_t = vertices[used], inverse.reshape(-1, 3)
    updated_keys, updated_ids = boundary_index(updated_v, updated_t)
    unresolved = [edges - updated_keys for edges in missing]
    return updated_v, updated_t, unresolved, updated_ids, np.flatnonzero(retry[~drop])


def stream_join(inputs, active, local_v, local_t, work):
    """Compact immutable interiors and append the patch using exact XY IDs.

    Input interiors occupy disjoint ownership boxes. Only patch coordinates
    can be shared across an interior/patch interface; their lookup is seam-sized.
    No global coordinate dictionary or global advanced-indexing copy is built.
    """
    mappings, counts = [], []
    for index, ((vertices, triangles, _), selection) in enumerate(zip(inputs, active)):
        ids = disk_array(work/f'ids-{index}.npy', (len(vertices),), np.int64)
        ids[:] = -1
        for lo, hi in batches(len(triangles)):
            ids[triangles[lo:hi][~selection[lo:hi]]] = 0
        count = 0
        for lo, hi in batches(len(ids)):
            used = np.flatnonzero(ids[lo:hi] >= 0)
            ids[lo+used] = sum(counts) + np.arange(count, count+len(used))
            count += len(used)
        counts.append(count)
        mappings.append(ids)
    lookup = {tuple(v[:2]): i for i, v in enumerate(local_v)}
    if len(lookup) != len(local_v):
        raise remesh.PatchError('Local merge contains coincident vertex IDs')
    local_ids = np.full(len(local_v), -1, dtype=np.int64)
    for (vertices, triangles, neighbours), selection, ids in zip(inputs, active, mappings):
        for _, endpoints, pair, _ in selected_boundary(vertices, triangles, neighbours, selection):
            for old, key in zip(endpoints, pair):
                local = lookup.get(key)
                if local is not None:
                    if local_ids[local] >= 0 and local_ids[local] != ids[old]:
                        raise remesh.PatchError('Immutable interiors share a vertex unexpectedly')
                    local_ids[local] = ids[old]
    added = np.flatnonzero(local_ids < 0)
    local_ids[added] = sum(counts) + np.arange(len(added))
    out_v = disk_array(work/'verts.npy', (sum(counts)+len(added), 3))
    nfaces = sum(sum(int(np.count_nonzero(~selection[lo:hi])) for lo, hi in batches(len(selection)))
                 for selection in active) + len(local_t)
    out_t = disk_array(work/'tris.npy', (nfaces, 3), np.int64)
    offset = 0
    for (vertices, triangles, _), selection, ids in zip(inputs, active, mappings):
        for lo, hi in batches(len(vertices)):
            used = lo+np.flatnonzero(ids[lo:hi] >= 0)
            out_v[ids[used]] = vertices[used]
        for lo, hi in batches(len(triangles)):
            block = ids[triangles[lo:hi][~selection[lo:hi]]]
            out_t[offset:offset+len(block)] = block
            offset += len(block)
    out_v[local_ids[added]] = local_v[added]
    for lo, hi in batches(len(local_t)):
        out_t[offset+lo:offset+hi] = local_ids[local_t[lo:hi]]
    out_v.flush(); out_t.flush()
    return out_v, out_t


def _write_cut_example(prefix, edge, parent, vertices, triangles, constraints):
    """Keep the local geometry of an exhausted cut without copying either tile.

    Include faces whose envelopes intersect the disk edge even if the seam
    omitted an endpoint and a matching subdivision chain cannot be identified.
    """
    local = []
    a, b = np.asarray(edge)
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    segment = LineString(edge)
    epsilon = remesh._coordinate_tolerance(segment)
    for start, end in batches(len(triangles)):
        xy = vertices[triangles[start:end], :2]
        nearby = ((xy.min(axis=1) <= hi+epsilon) & (xy.max(axis=1) >= lo-epsilon)).all(axis=1)
        for tri in xy[nearby]:
            if Polygon(tri).distance(segment) <= epsilon:
                local.append(tri.tolist())
    envelope = union_all([Polygon(parent), *(Polygon(tri) for tri in local)]).envelope
    lines = [list(part.coords) for curve in constraints for part in remesh.line_parts(curve.intersection(envelope))]
    path = Path(str(prefix) + '.unmatched-edge.json')
    write_marker(path, dict(disk_edge=edge, disk_parent=parent.tolist(), seam_triangles=local,
                            constraint_lines=lines))
    return path


def _save_retry(path, active, pending, cached, roundoff, state):
    """Commit the *next* outer attempt after selecting its new disk collar.

    Parent meshes and neighbour tables can be reopened by ID. Only the two
    selection masks, the accepted local mesh and its collar mappings need to
    survive the PBS allocation. The previous generation stays usable up to rename.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {'state': json.dumps(state, default=lambda value: value.item())}
    for side in (0, 1):
        arrays[f'active{side}'] = np.packbits(active[side])
        arrays[f'roundoff{side}'] = np.array(sorted(roundoff[side]), dtype=np.int64)
        if cached is not None:
            arrays[f'added{side}'], arrays[f'removed{side}'] = pending[side]
            arrays[f'shared{side}'] = np.array(list(cached[2][side].items()), dtype=np.int64).reshape(-1, 2)
    if cached is not None:
        arrays['vertices'], arrays['triangles'], arrays['release'] = cached[0], cached[1], cached[3]
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        tmp = Path(stream.name)
        try:
            np.savez(stream, **arrays)
            stream.close()
            os.replace(tmp, path)
        finally:
            tmp.unlink(missing_ok=True)


def _load_retry(path, active):
    """Read and CRC-check the entire snapshot before mutating either input mask."""
    try:
        with np.load(path, allow_pickle=False) as data:
            state = json.loads(data['state'].item())
            masks = [np.unpackbits(data[f'active{i}'], count=len(selection))
                     for i, selection in enumerate(active)]
            roundoff = [set(data[f'roundoff{i}']) for i in (0, 1)]
            cached = pending = None
            if 'vertices' in data:
                shared = [dict(data[f'shared{i}']) for i in (0, 1)]
                # Older checkpoints retained the complete local candidate.
                release = data['release'] if 'release' in data else np.empty(0, dtype=np.int64)
                cached = data['vertices'], data['triangles'], shared, release
                pending = [(data[f'added{i}'], data[f'removed{i}']) for i in (0, 1)]
    except FileNotFoundError:
        return None
    except (ValueError, KeyError, EOFError, zipfile.BadZipFile, zlib.error):
        return None
    for selection, mask in zip(active, masks):
        selection[:] = mask
    return cached, pending, roundoff, state


# Remesh the seam while the two interiors stay on disk. Compare the solved boundary against both
# disk cuts before joining them. On a mismatch, carry the solved faces into a wider collar
# within the active-face budget.
def merge_disk_cavity(path_a, bounds_a, path_b, bounds_b, domain, constraints, options, output):
    started = time.perf_counter()
    # Keep failed native reproductions beside the final mesh even when staging its archive.
    directory = str(Path(options.get('out_prefix') or output).parent)
    width = options.get('remesh_band_width')
    width = float(2*math.sqrt(options['max_area']) if width is None else width)
    limit = int(options.get('remesh_max_attempts', 12))
    max_active = int(options.get('remesh_max_active_faces', 2000000))
    if limit < 1 or max_active < 1 or not np.isfinite(width) or width <= 0:
        raise remesh.PatchError('Seam width, attempt count and active-face limit must be positive')
    axis, seam = remesh._shared_interface(bounds_a, bounds_b)
    bounds = [min(bounds_a[0], bounds_b[0]), min(bounds_a[1], bounds_b[1]),
              max(bounds_a[2], bounds_b[2]), max(bounds_a[3], bounds_b[3])]
    target = remesh._clipped_domain(domain, bounds, 0, width/2, options.get('out_prefix'),
                                    options.get('mesher_debug', False))
    with tempfile.TemporaryDirectory(prefix='seam-store-', dir=directory) as name, ExitStack() as stack:
        work = Path(name)
        inputs, active = [], []
        input_seconds = 0.0
        input_neighbour_seconds = 0.0
        for index, (path, core) in enumerate(((path_a, bounds_a), (path_b, bounds_b))):
            report_phase(options, f'opening input {index+1}/2')
            read_started = time.perf_counter()
            v, t = stack.enter_context(open_mesh_arrays(path, work))
            input_seconds += time.perf_counter() - read_started
            report_phase(options, f'building neighbours for input {index+1}/2 ({len(t):,} faces)')
            neighbour_started = time.perf_counter()
            neighbours = (build_neighbours(t, len(v), v, chunk_size=32768,
                output_path=work/f'neighbours-{index}.npy', temporary_directory=work)
                if len(t) else np.empty((0, 3), dtype=np.int64))
            input_neighbour_seconds += time.perf_counter() - neighbour_started
            # Active faces belong to the local solve. This includes collar
            # triangles we hope to keep, so the disk cut can move outward
            # without forcing every newly loaded face to be remeshed.
            report_phase(options, f'selecting seam faces in input {index+1}/2')
            selection = disk_array(work/f'active-{index}.npy', (len(t),), bool)
            for lo, hi in batches(len(t)):
                selection[lo:hi] = ~remesh._retain_mask(v, t[lo:hi], core, axis, seam, width)
            # A four-ring collar gives local parent repairs room before touching disk interiors.
            seeds = []
            for lo, hi in batches(len(t)):
                adjacent = neighbours[lo:hi]
                touches_active = ((adjacent >= 0) & selection[adjacent]).any(axis=1)
                seeds.extend(lo + np.flatnonzero(~selection[lo:hi] & touches_active))
            grow_active(selection, neighbours, seeds, 3)
            inputs.append((v, t, neighbours)); active.append(selection)
        reports = []
        cached = None
        pending = None
        total = sum(len(t) for _, t, _ in inputs)
        peak_active = 0
        reused = 0
        roundoff_parents = [set(), set()]

        # A roundoff-width strip can be bounded by a face outside the active mesh. Find that
        # face on the disk cut and move it into the solve that can change its edge.
        def release_disk_parents(chords):
            count = 0
            for values, selection, released in zip(inputs, active, roundoff_parents):
                boundary = selected_boundary(*values, selection)
                seeds = remesh._roundoff_owners(chords, ((face, xy) for face, _, xy, _ in boundary))
                selection[seeds] = True
                released.update(seeds)
                count += len(seeds)
            if count:
                debug_log(options, f'Remeshing {count} disk faces along a roundoff-width boundary')
            return count

        previous_unmatched = None
        collar_rings = 4
        attempt = 0
        retry_path = options.get('_seam_retry_path')
        saved = _load_retry(retry_path, active) if retry_path is not None else None
        if saved is not None:
            cached, pending, roundoff_parents, state = saved
            attempt, collar_rings, previous_unmatched, reused, peak_active, reports = state
            print(f'Resuming seam attempt {attempt+1}/{limit} [{Path(output).stem}]', flush=True)
            del saved
        while attempt < limit:
            report_phase(options, f'finding disk interfaces for seam attempt {attempt+1}/{limit}')
            # The mask counts original faces for disk ownership. After a retry
            # their replacements may contain many more triangles, so charge
            # the cached triangle count against the next solve's memory budget.
            selected_count = sum(int(np.count_nonzero(selection)) for selection in active)
            active_count = (selected_count if cached is None
                            else len(cached[1]) + sum(len(added) for added, _ in pending))
            coverage, cuts, cut_vertices = [], [], []
            for side, (values, selection) in enumerate(zip(inputs, active)):
                polygon, interfaces, endpoints = inactive_coverage(*values, selection)
                # The last solve can split a previously matching edge on either
                # tile. Close the remaining cut-connected interiors that fit
                # after reserving the cached candidate and its pending collars.
                if attempt > 0 and attempt+1 == limit and interfaces and active_count < max_active:
                    seeds = np.unique(list(interfaces.values()))
                    added = grow_active(selection, values[2], seeds, len(selection),
                                        max_faces=max_active-active_count)
                    if added is not None and len(added):
                        if cached is not None:
                            collar, removed = pending[side]
                            pending[side] = (np.concatenate((collar, added)), removed)
                        selected_count += len(added)
                        active_count += len(added)
                        debug_log(options, f'Final seam recovery [{Path(output).stem}]: '
                                  f'loaded {len(added)} faces touching the disk interface')
                        polygon, interfaces, endpoints = inactive_coverage(*values, selection)
                coverage.append(polygon)
                cuts.append(interfaces)
                cut_vertices.append(endpoints)
            outside = remesh.polygon_domain(union_all(coverage))
            if outside.difference(target).area > remesh.area_tolerance(target):
                raise remesh.PatchError('Immutable interiors extend beyond the merge domain')
            local_domain = remesh.polygon_domain(target.difference(outside))
            peak_active = max(peak_active, active_count)
            if active_count > max_active:
                raise remesh.PatchError(f'Seam active-face budget exceeded: {active_count} > {max_active}; '
                                       'increase mpi_remesh_max_active_faces only if memory permits')
            if not outside.is_empty:
                chords = [*remesh._roundoff_chords(local_domain),
                          *remesh._unbisectable_edges(remesh.line_parts(local_domain.boundary))]
                if chords and release_disk_parents(chords):
                    cached = None
                    continue
            if cached is None:
                av, at, a_ids = compact_region(*inputs[0][:2], active[0])
                bv, bt, b_ids = compact_region(*inputs[1][:2], active[1])
                initial_keep = None
                if any(roundoff_parents):
                    keep_a = remesh._retain_mask(av, at, bounds_a, axis, seam, width)
                    keep_b = remesh._retain_mask(bv, bt, bounds_b, axis, seam, width)
                    initial_keep = np.concatenate((keep_a & ~np.isin(a_ids, list(roundoff_parents[0])),
                                                   keep_b & ~np.isin(b_ids, list(roundoff_parents[1]))))
            else:
                reused += len(cached[1]) - len(cached[3])
                av, at, initial_keep = extend_cached_region(cached, inputs, pending)
                # The expanded seam is already indexed as a single working mesh.
                bv, bt = av[:0], at[:0]
            debug_log(options, f'Seam region [{Path(output).stem}]: active={active_count}, '
                      f'interiors_on_disk={total-selected_count}')
            report_phase(options, f'solving seam, attempt {attempt+1}/{limit} ({active_count:,} active faces)')
            # One disk attempt can spend hours in its own series of cavity
            # solves. Keep the outer counter visible on every inner heartbeat.
            seam_phase = f'seam attempt {attempt+1}/{limit}'
            if previous_unmatched is not None:
                seam_phase += f' (unmatched disk edges at last check = {previous_unmatched:,})'
            try:
                lv, lt, report = remesh.merge_cavity(av, at, bounds_a, bv, bt, bounds_b,
                    domain, constraints, options, directory, prepared_domain=local_domain, initial_keep=initial_keep,
                    disk_coverage=outside, phase_prefix=seam_phase + ' | ')
                report_phase(options, f'checking disk interfaces for seam attempt {attempt+1}/{limit}')
                # Splitting a disk-facing edge leaves the covered area unchanged.
                # Check its endpoints too, or the split vertex becomes a hanging
                # node when the two meshes are joined.
                if len(lt):
                    # Both checks inspect the same accepted triangles. Share their
                    # adjacency here, before an interface repair can change face IDs.
                    local_neighbours = build_neighbours(lt, len(lv), lv)
                    local_coverage = remesh.mesh_coverage(lv, lt, local_neighbours)
                    if local_coverage.symmetric_difference(local_domain).area > remesh.area_tolerance(local_domain):
                        raise remesh.PatchError('Seam region has inconsistent target coverage')
                    keys, boundary_ids = boundary_index(lv, lt, local_neighbours)
                    del local_neighbours
                elif not local_domain.is_empty:
                    raise remesh.PatchError('Seam region is unexpectedly empty')
                else:
                    keys, boundary_ids = set(), {}
                missing = [interfaces.keys() - keys for interfaces in cuts]
                reports.append(report)
                report['unmatched_disk_edges'] = sum(map(len, missing))
                affected = [sorted({interfaces[edge] for edge in edges}) for interfaces, edges in zip(cuts, missing)]
                parent_count = sum(map(len, affected))
                if parent_count and limit > 1 and active_count + parent_count <= max_active:
                    peak_active = max(peak_active, active_count + parent_count)
                    repaired = subdivide_disk_interface(inputs, active, affected, missing, lv, lt,
                                                        boundary_ids, keys, options, directory)
                    if repaired is not None:
                        lv, lt, accepted, missing, cut_vertices, boundary_ids, keys = repaired
                        for selection, faces in zip(active, accepted):
                            selection[faces] = True
                        accepted_count = sum(map(len, accepted))
                        selected_count += accepted_count
                        affected = [sorted({interfaces[e] for e in edges}) for interfaces, edges in zip(cuts, missing)]
                        parent_count = sum(map(len, affected))
                        report['disk_parents_subdivided'] = accepted_count
                        debug_log(options, f'Subdivided {accepted_count} disk parents [{Path(output).stem}]: '
                                  f'{parent_count} parents still require expansion')
                release = np.empty(0, dtype=np.int64)
                remaining = sum(map(len, missing))
                stalled = previous_unmatched is not None and remaining >= previous_unmatched
                if remaining and limit > 1 and (stalled or attempt+1 == limit):
                    report_phase(options, f'repairing {remaining} disk edges from the seam side')
                    lv, lt, missing, boundary_ids, release = repair_seam_boundary(
                        lv, lt, keys, boundary_ids, missing, cut_vertices, constraints, options, directory)
                    report['disk_edges_restored_from_seam'] = remaining - sum(map(len, missing))
                    report['seam_faces_released_for_retry'] = len(release)
                    affected = [sorted({interfaces[e] for e in edges}) for interfaces, edges in zip(cuts, missing)]
                    parent_count = sum(map(len, affected))
                    debug_log(options, f"Restored {report['disk_edges_restored_from_seam']} disk edges "
                              f"from the seam side; {parent_count} disk parents still unmatched")
                report['unresolved_disk_edges'] = sum(map(len, missing))
                if any(affected):
                    shared = [{i: boundary_ids[xy] for i, xy in endpoints.items() if xy in boundary_ids}
                              for endpoints in cut_vertices]
                    cached = lv, lt, shared, release
            except remesh.PatchBoundaryRoundoff as exc:
                # Preflight can release another disk parent before CGAL runs.
                # Each successful release shrinks the immutable interior, so
                # this path can make progress even on the last allowed solve.
                if not exc.preflight_only:
                    attempt += 1
                if attempt == limit or not release_disk_parents(exc.chords):
                    raise
                cached = None
                continue
            except remesh.PhysicalConstraintAngle:
                raise
            except remesh.PatchError as exc:
                cause = exc
                while cause is not None:
                    if (isinstance(cause, remesh.BackendTerminated) or 'Patch vertex budget exceeded' in str(cause)
                            or 'Seam active-face budget exceeded' in str(cause)):
                        raise
                    cause = cause.__cause__
                # A strict quality failure at an artificial outer cut can require
                # more neighbours. Never relax the criteria to keep an interior.
                affected = [sorted(set(interfaces.values())) for interfaces in cuts]
                if not any(affected) or attempt+1 == limit:
                    raise
                cached = None
                debug_log(options, f'Expanding seam region after rejected local repair: {exc}')
            else:
                if not any(affected):
                    report_phase(options, 'joining the seam to the disk interiors')
                    out_v, out_t = stream_join(inputs, active, lv, lt, work)
                    output_neighbour_seconds = 0.0
                    if len(out_t):
                        report_phase(options, f'building output neighbours ({len(out_t):,} faces)')
                        neighbour_started = time.perf_counter()
                        neighbours = build_neighbours(out_t, len(out_v), out_v, 32768,
                            work/'final-neighbours.npy', work)
                        output_neighbour_seconds = time.perf_counter() - neighbour_started
                        report_phase(options, 'checking the merged boundary')
                        check_boundary_geometry(out_v, out_t, neighbours, work)
                    # Copy NPY members directly into the backward-compatible NPZ archive.
                    report_phase(options, 'writing the merged archive')
                    write_started = time.perf_counter()
                    temporary = work/'complete.npz'
                    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_STORED, allowZip64=True) as archive:
                        archive.write(work/'verts.npy', 'verts.npy')
                        archive.write(work/'tris.npy', 'tris.npy')
                    os.replace(temporary, output)
                    report.update(seam_attempts=attempt+1, seam_attempt_reports=reports[:-1],
                        peak_loaded_faces=peak_active, seam_work_reused=reused,
                        roundoff_disk_parents=sum(map(len, roundoff_parents)),
                        boundary_halo_area=target.area-domain.intersection(box(*bounds)).area,
                        disk_interior_faces=len(out_t)-len(lt), output_faces=len(out_t),
                        input_seconds=input_seconds, input_neighbour_seconds=input_neighbour_seconds,
                        output_neighbour_seconds=output_neighbour_seconds, output_npz_seconds=time.perf_counter()-write_started,
                        elapsed_seconds=time.perf_counter()-started)
                    return report
            if attempt+1 == limit:
                side = next(i for i, edges in enumerate(missing) if edges)
                edge = min(missing[side])
                example = tuple(tuple(map(float, p)) for p in edge)
                remaining = [r['unresolved_disk_edges'] for r in reports if 'unresolved_disk_edges' in r]
                v, t, _ = inputs[side]
                witness = _write_cut_example(options.get('out_prefix') or Path(directory)/Path(output).stem,
                                             example, v[t[cuts[side][edge]], :2], lv, lt, constraints)
                raise remesh.PatchError(
                    f'Seam-region expansion exhausted [{Path(output).stem}] after {attempt+1} attempts: '
                    f'{sum(map(len, missing))} unmatched disk edges on {parent_count} parents, '
                    f'active={active_count}, active_face_limit={max_active}, '
                    f'interiors_on_disk={total-selected_count}, unresolved_by_solve={remaining}, example={example}, '
                    f'owner={cuts[side][edge]} in {(path_a, path_b)[side]}, inputs=({path_a}, {path_b}), '
                    f'local geometry={witness}')
            pending = []
            for selection, (_, _, neighbours), seeds in zip(active, inputs, affected):
                # Remove the failed interface parents and their immediate ring;
                # load three further rings as retained repair candidates.
                if seeds:
                    adjacent = neighbours[seeds].ravel()
                    removed = np.unique(np.concatenate((seeds, adjacent[adjacent >= 0])))
                    removed = removed[~selection[removed]]
                    added = grow_active(selection, neighbours, seeds, 4)
                else:
                    removed = added = np.empty(0, dtype=np.int64)
                pending.append((added, removed))
            if cached is not None:
                # Reserve both ordinary collars first. On the last expansion,
                # remove the artificial cut where its whole component fits.
                capacity = max_active - len(lt) - sum(len(added) for added, _ in pending)
                close_components = attempt+2 == limit and close_disk_components(inputs, active, pending, capacity)
                # A repair can split the next disk edge and move the same problem
                # outward one ring at a time. Back off by doubling the retained
                # collar while the count stalls. Once it falls, local growth
                # saves memory. Only the failed parents and their immediate
                # neighbours are selected for removal here.
                unmatched = report['unresolved_disk_edges']
                if previous_unmatched is not None and unmatched >= previous_unmatched:
                    collar_rings = min(2 * collar_rings, max_active)
                else:
                    collar_rings = 4
                previous_unmatched = unmatched
                if collar_rings > 4 and not close_components:
                    widen_disk_collars(inputs, active, pending, capacity, collar_rings - 4)
                count = sum(len(added) for added, _ in pending)
                debug_log(options, f'Expanding disk interface [{Path(output).stem}]: '
                          f'unmatched_edges={unmatched}, loading={count}, '
                          f'closing_components={close_components}')
                report['active_faces'] = active_count
                report['disk_expansion_faces'] = count
                report['closing_disk_components'] = close_components
            else:
                previous_unmatched = None
            attempt += 1
            if retry_path is not None:
                _save_retry(retry_path, active, pending, cached, roundoff_parents,
                            [attempt, collar_rings, previous_unmatched, reused, peak_active, reports])
        raise remesh.PatchError('No conforming seam-region replacement')
