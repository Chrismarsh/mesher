"""Triangle connectivity shared by MPI export and mesh inspection tools.

Neighbour slot i is opposite vertex i. Edges are matched by vertex ID, never
by coordinate proximity. Large inputs use sorted temporary runs so edge
matching does not require a Python dictionary containing the whole mesh.
"""

from contextlib import ExitStack
import os
import tempfile

import numpy as np


class MeshTopologyError(ValueError):
    """The mesh cannot represent a valid triangular neighbour graph."""


EDGE_DTYPE = np.dtype([('a', '<i8'), ('b', '<i8'), ('owner', '<i8')])
DEFAULT_CHUNK_SIZE = 1_000_000


# Connectivity is indexed directly into arrays. Check the shape and range before converting to
# signed IDs so an oversized unsigned index cannot wrap to a boundary sentinel.
def integer_rows(values, name, width=3):
    array = np.asarray(values)
    if array.ndim != 2 or array.shape[1] != width:
        raise MeshTopologyError(f'{name} must have shape (N, {width}); got {array.shape}')
    if array.dtype.kind not in 'iu':
        raise MeshTopologyError(f'{name} must contain integer indices')
    if array.dtype.kind == 'u' and array.size and array.max() > np.iinfo(np.int64).max:
        raise MeshTopologyError(f'{name} contains indices exceeding signed 64-bit range')
    return np.asarray(array, dtype=np.int64)


# Sort edges by their unordered endpoint IDs so owners from different faces meet. Packing face *
# 3 + slot lets the match write directly into the neighbour table.
def _edge_records(triangles, first_face):
    edges = np.empty(3 * len(triangles), dtype=EDGE_DTYPE)
    for slot, (a, b) in enumerate(((1, 2), (2, 0), (0, 1))):
        edges['a'][slot::3] = np.minimum(triangles[:, a], triangles[:, b])
        edges['b'][slot::3] = np.maximum(triangles[:, a], triangles[:, b])
    edges['owner'] = np.arange(3 * first_face, 3 * (first_face + len(triangles)))
    # Owners already increase in input order. The stable lexsort also leaves equal
    # endpoint pairs in owner order, as required by the on-disk merge.
    return edges[np.lexsort((edges['b'], edges['a']))]


# One owner leaves a boundary, two owners fill reciprocal slots, and a third makes adjacency
# ambiguous. Equal endpoint pairs are consecutive after the sort.
def _match_sorted_edges(edges, neighbours):
    same = (edges['a'][1:] == edges['a'][:-1]) & (edges['b'][1:] == edges['b'][:-1])
    triples = np.flatnonzero(same[:-1] & same[1:])
    if len(triples):
        edge = edges[triples[0]]
        raise MeshTopologyError(f'Non-manifold edge ({edge["a"]}, {edge["b"]}): more than two triangles')
    pairs = np.flatnonzero(same)
    left = edges['owner'][pairs]
    right = edges['owner'][pairs + 1]
    neighbours.flat[left] = right // 3
    neighbours.flat[right] = left // 3


def _merged_edge_blocks(paths, block_size):
    """Yield sorted blocks from sorted scratch files, with a bounded working set.

    Use the complete (a, b, owner) key. The owner is unique even for a
    non-manifold edge, so a dense vertex fan or corrupt input cannot push
    thousands of equal keys past the per-run read limit.
    """
    with ExitStack() as mappings:
        runs = []
        for path in paths:
            run = np.memmap(path, dtype=EDGE_DTYPE, mode='r')
            mappings.callback(run._mmap.close)
            runs.append(run)

        while len(runs) > 1:
            # Each run can supply at most its share of this block. Up to the
            # smallest end key, all records are present in those prefixes.
            share = block_size // len(runs)
            pivot = min((run[share - 1] for run in runs if len(run) > share),
                        key=tuple, default=None)
            if pivot is None:
                counts = [len(run) for run in runs]
            else:
                counts = [run[:share].searchsorted(pivot, side='right') for run in runs]
            parts = [run[:n] for run, n in zip(runs, counts) if n]
            runs = [run[n:] for run, n in zip(runs, counts) if n < len(run)]
            if len(parts) == 1:
                yield parts[0]
            else:
                records = np.concatenate(parts)
                yield records[np.lexsort((records['owner'], records['b'], records['a']))]

        if runs:
            for start in range(0, len(runs[0]), block_size):
                yield runs[0][start:start+block_size]


def _match_edge_runs(paths, neighbours, directory, block_size):
    """Pair owners across runs without visiting every edge in a Python loop."""
    # Limit open files even when a late MPI stage has hundreds of runs. Each
    # remaining file gets at least one record in a merge block.
    fan_in = min(64, block_size)
    generation = 0
    while len(paths) > fan_in:
        next_paths = []
        for start in range(0, len(paths), fan_in):
            group = paths[start:start+fan_in]
            if len(group) == 1:
                next_paths.extend(group)
                continue
            path = os.path.join(directory, f'merged-{generation}-{start}')
            with open(path, 'wb') as output:
                for block in _merged_edge_blocks(group, block_size):
                    block.tofile(output)
            for old in group:
                os.unlink(old)
            next_paths.append(path)
        paths = next_paths
        generation += 1

    tail = np.empty(0, dtype=EDGE_DTYPE)
    for block in _merged_edge_blocks(paths, block_size):
        # A pair or a non-manifold triple can straddle a block boundary.
        # Rechecking the last two records also covers successive one-record blocks.
        edges = np.concatenate((tail, block))
        _match_sorted_edges(edges, neighbours)
        tail = edges[-2:].copy()


def build_neighbours(elements, nvertex, vertices=None, chunk_size=DEFAULT_CHUNK_SIZE,
                     output_path=None, temporary_directory=None):
    """Return an int64 (N, 3) table, rejecting invalid triangle topology.

    chunk_size bounds the number of triangles sorted in memory at once.
    An output_path keeps the neighbour table on disk. Edge sorting and geometry
    checks use bounded batches; temporary_directory selects the scratch volume.
    """
    triangles = integer_rows(elements, 'elem')
    if not len(triangles):
        raise MeshTopologyError('Mesh has no triangles')
    if not isinstance(chunk_size, (int, np.integer)) or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    chunk_size = int(chunk_size)  # Three edge records per face can overflow a narrow NumPy scalar.
    if triangles.min() < 0 or triangles.max() >= nvertex:
        raise MeshTopologyError(f'elem contains a vertex index outside [0, {nvertex})')
    positions = None
    if vertices is not None:
        positions = np.asarray(vertices, dtype=float)
        if positions.ndim != 2 or positions.shape[0] != nvertex or positions.shape[1] not in (2, 3):
            raise MeshTopologyError('vertex must have shape (nvertex, 2) or (nvertex, 3)')
        for start in range(0, len(positions), chunk_size):
            if not np.isfinite(positions[start:start + chunk_size]).all():
                raise MeshTopologyError('vertex contains non-finite coordinates or elevations')
    for start in range(0, len(triangles), chunk_size):
        block = triangles[start:start + chunk_size]
        repeated = (block[:, 0] == block[:, 1]) | (block[:, 1] == block[:, 2]) | (block[:, 0] == block[:, 2])
        if repeated.any():
            raise MeshTopologyError(f'Triangle {start + np.flatnonzero(repeated)[0]} repeats a vertex')
        if positions is not None:
            u = positions[block[:, 1], :2] - positions[block[:, 0], :2]
            v = positions[block[:, 2], :2] - positions[block[:, 0], :2]
            p, q = u[:, 0] * v[:, 1], u[:, 1] * v[:, 0]
            degenerate = np.abs(p - q) <= 16 * np.finfo(float).eps * (np.abs(p) + np.abs(q))
            if degenerate.any():
                raise MeshTopologyError(f'Triangle {start + np.flatnonzero(degenerate)[0]} has zero or numerically degenerate area')

    neighbours = (np.empty(triangles.shape, dtype=np.int64) if output_path is None else
                  np.lib.format.open_memmap(output_path, mode='w+', dtype=np.int64, shape=triangles.shape))
    neighbours[:] = -1
    if len(triangles) <= chunk_size:
        _match_sorted_edges(_edge_records(triangles, 0), neighbours)
    else:
        with tempfile.TemporaryDirectory(prefix='mesher-edges-', dir=temporary_directory) as directory:
            paths = []
            for start in range(0, len(triangles), chunk_size):
                path = os.path.join(directory, str(start))
                _edge_records(triangles[start:start + chunk_size], start).tofile(path)
                paths.append(path)
            _match_edge_runs(paths, neighbours, directory, 3 * chunk_size)

    for start in range(0, len(neighbours), chunk_size):
        block = neighbours[start:start + chunk_size]
        duplicates = ((block[:, 0] >= 0) & ((block[:, 0] == block[:, 1]) | (block[:, 0] == block[:, 2]))) | ((block[:, 1] >= 0) & (block[:, 1] == block[:, 2]))
        if duplicates.any():
            raise MeshTopologyError(f'Duplicate triangles involving face {start + np.flatnonzero(duplicates)[0]}')
    return neighbours


def validate_neighbours(values, expected=None):
    """Validate a neighbour graph and optionally its exact edge-slot mapping."""
    if len(values) == 0:
        raise MeshTopologyError('Neighbour table is empty. Audit geometry before rebuilding neighbours or partitioning')
    neighbours = integer_rows(values, 'neigh')
    n = len(neighbours)
    if expected is not None and neighbours.shape != expected.shape:
        raise MeshTopologyError(f'Expected {len(expected)} neighbour rows; got {n}')
    if neighbours.min() < -1 or neighbours.max() >= n:
        raise MeshTopologyError(f'neigh indices must be -1 or in [0, {n})')
    for start in range(0, n, 100_000):
        block = neighbours[start:start + 100_000]
        ids = np.arange(start, start + len(block))
        if (block == ids[:, None]).any():
            raise MeshTopologyError('neigh contains a self-neighbour')
        for slot in range(3):
            valid = block[:, slot] >= 0
            targets = block[valid, slot]
            if not np.all(np.any(neighbours[targets] == ids[valid, None], axis=1)):
                raise MeshTopologyError('neigh is not reciprocal')
        repeated = ((block[:, 0] >= 0) & ((block[:, 0] == block[:, 1]) | (block[:, 0] == block[:, 2]))) | ((block[:, 1] >= 0) & (block[:, 1] == block[:, 2]))
        if repeated.any():
            raise MeshTopologyError('neigh contains a repeated neighbour')
        if expected is not None and not np.array_equal(block, expected[start:start + len(block)]):
            raise MeshTopologyError('neigh disagrees with triangle edges or opposite-vertex slots')
    return neighbours


def boundary_subdivisions(vertices, triangles, neighbours, tolerance):
    """Return coincidence count and an iterator over boundary-edge subdivisions.

    Each item is (face, slot, ordered_vertex_ids). Only boundary vertices strictly
    inside an unmatched edge, within tolerance in XY, are considered. This is a
    geometric candidate search, not proof that the edge is an MPI seam.
    """
    coincident, subdivisions = boundary_node_conflicts(vertices, triangles, neighbours, tolerance)
    return len(coincident), subdivisions


def boundary_node_conflicts(vertices, triangles, neighbours, tolerance):
    """Return coincident boundary vertex IDs and the same edge subdivision iterator."""
    from scipy.spatial import cKDTree

    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError('boundary tolerance must be finite and non-negative')
    faces, slots = np.nonzero(neighbours == -1)
    a = triangles[faces, (slots + 1) % 3]
    b = triangles[faces, (slots + 2) % 3]
    ids = np.unique(np.concatenate((a, b)))
    if len(ids) < 2:
        return ids[:0], iter(())
    xy = vertices[ids, :2]
    tree = cKDTree(xy)
    distances, _ = tree.query(xy, k=2)
    coincident = ids[distances[:, 1] <= tolerance]

    # Project candidate vertices onto each open segment. Returning interior IDs in distance
    # order gives the stitcher a chain from the first endpoint to the second.
    def subdivisions():
        for face, slot, va, vb in zip(faces, slots, a, b):
            p, q = vertices[va, :2], vertices[vb, :2]
            delta = q - p
            length = np.linalg.norm(delta)
            candidates = tree.query_ball_point((p + q) / 2, length / 2 + tolerance)
            inside = []
            for candidate in candidates:
                vid = int(ids[candidate])
                if vid in (va, vb):
                    continue
                offset = xy[candidate] - p
                along = np.dot(offset, delta) / length
                distance = abs(delta[0] * offset[1] - delta[1] * offset[0]) / length
                if tolerance < along < length - tolerance and distance <= tolerance:
                    inside.append((along, vid))
            if inside:
                yield int(face), int(slot), [vid for _, vid in sorted(inside)]

    return coincident, subdivisions()


def _boundary_defects(vertices, triangles, neighbours, tolerance):
    """Find coincident boundary vertices and T-junctions; do not weld them."""
    coincident, subdivisions = boundary_subdivisions(vertices, triangles, neighbours, tolerance)
    hanging = 0
    examples = []
    for face, slot, nodes in subdivisions:
        va, vb = triangles[face, (slot + 1) % 3], triangles[face, (slot + 2) % 3]
        hanging += len(nodes)
        for vid in nodes[:max(0, 10 - len(examples))]:
            examples.append({'vertex': vid, 'edge': [int(va), int(vb)], 'face': face,
                             'vertex_xy': vertices[vid, :2].tolist(),
                             'edge_xy': vertices[[va, vb], :2].tolist()})
    return {'coincident_boundary_vertices': coincident, 'hanging_nodes': hanging, 'hanging_node_examples': examples}


def audit_mesh(document, chunk_size=DEFAULT_CHUNK_SIZE, boundary_tolerance=1e-8, check_boundaries=True):
    """Return (report, reconstructed_neighbours); never modify document.

    This checks indexed topology and boundary T-junctions, not arbitrary polygon
    overlaps or coverage of the intended model domain. Disconnected components
    are reported because they can be intentional.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    if not np.isfinite(boundary_tolerance) or boundary_tolerance < 0:
        raise ValueError('boundary_tolerance must be finite and non-negative')
    report = {'errors': [], 'neighbour_errors': [], 'warnings': [], 'neighbour_rebuild_allowed': False,
              'boundary_geometry_checked': check_boundaries, 'domain_coverage_checked': False}
    try:
        mesh = document['mesh']
        vertices = np.asarray(mesh['vertex'], dtype=float)
        triangles = integer_rows(mesh['elem'], 'elem')
        for name, actual in (('nvertex', len(vertices)), ('nelem', len(triangles))):
            if not isinstance(mesh.get(name), int) or isinstance(mesh.get(name), bool) or mesh[name] != actual:
                raise MeshTopologyError(f'{name} must equal the actual row count ({actual})')
        if vertices.ndim != 2 or vertices.shape[1] != 3:
            raise MeshTopologyError('vertex must have shape (nvertex, 3)')
        neighbours = build_neighbours(triangles, len(vertices), vertices, chunk_size)
    except (KeyError, TypeError, ValueError) as exc:
        report['errors'].append(str(exc))
        report['valid'] = False
        return report, None

    report.update(nvertex=len(vertices), nelem=len(triangles), boundary_edges=int(np.count_nonzero(neighbours == -1)))
    try:
        validate_neighbours(mesh.get('neigh', []), neighbours)
    except (TypeError, ValueError) as exc:
        report['neighbour_errors'].append(str(exc))

    indices = neighbours[neighbours >= 0]
    indptr = np.empty(len(neighbours) + 1, dtype=np.int64)
    indptr[0] = 0
    np.cumsum(np.sum(neighbours >= 0, axis=1), out=indptr[1:])
    graph = csr_matrix((np.ones(len(indices), dtype=np.uint8), indices, indptr), shape=(len(neighbours), len(neighbours)))
    count, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels)
    report['components'] = int(count)
    report['largest_component_sizes'] = np.sort(sizes)[-10:][::-1].tolist()
    if count > 1:
        report['warnings'].append(f'{count} connected components; verify these are intentional')
    if check_boundaries:
        defects = _boundary_defects(vertices, triangles, neighbours, boundary_tolerance)
        report.update(defects)
        if defects['coincident_boundary_vertices']:
            report['errors'].append('Coincident boundary vertices use different IDs; neighbour-only repair cannot join these seams')
        if defects['hanging_nodes']:
            report['errors'].append('Hanging boundary nodes detected; neighbour-only repair cannot fix these seams')
    report['valid'] = not (report['errors'] or report['neighbour_errors'])
    report['neighbour_rebuild_allowed'] = check_boundaries and not report['errors']
    if report['errors']:
        report['recovery'] = 'Geometry must be corrected before neighbour-only recovery; --rebuild-neighbours is blocked'
    elif report['neighbour_errors']:
        report['recovery'] = 'Run meshaudit --rebuild-neighbours -o NEW.mesh to rebuild neighbours without changing triangles'
    return report, neighbours
