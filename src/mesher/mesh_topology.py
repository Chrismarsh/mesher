"""Triangle connectivity shared by MPI export and mesh inspection tools.

Neighbour slot i is opposite vertex i. Edges are matched by vertex ID, never
by coordinate proximity. Large inputs use sorted temporary runs so edge
matching does not require a Python dictionary containing the whole mesh.
"""

import heapq
import os
import tempfile

import numpy as np


class MeshTopologyError(ValueError):
    """The mesh cannot represent a valid triangular neighbour graph."""


EDGE_DTYPE = np.dtype([('a', '<i8'), ('b', '<i8'), ('owner', '<i8')])
DEFAULT_CHUNK_SIZE = 1_000_000


def integer_rows(values, name, width=3):
    array = np.asarray(values)
    if array.ndim != 2 or array.shape[1] != width:
        raise MeshTopologyError(f'{name} must have shape (N, {width}); got {array.shape}')
    if array.dtype.kind not in 'iu':
        raise MeshTopologyError(f'{name} must contain integer indices')
    if array.dtype.kind == 'u' and array.size and array.max() > np.iinfo(np.int64).max:
        raise MeshTopologyError(f'{name} contains indices exceeding signed 64-bit range')
    return np.asarray(array, dtype=np.int64)


def _edge_records(triangles, first_face):
    edges = np.empty(3 * len(triangles), dtype=EDGE_DTYPE)
    for slot, (a, b) in enumerate(((1, 2), (2, 0), (0, 1))):
        edges['a'][slot::3] = np.minimum(triangles[:, a], triangles[:, b])
        edges['b'][slot::3] = np.maximum(triangles[:, a], triangles[:, b])
    edges['owner'] = np.arange(3 * first_face, 3 * (first_face + len(triangles)))
    edges.sort(order=['a', 'b', 'owner'])
    return edges


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


def _records(array):
    for record in array:
        yield int(record['a']), int(record['b']), int(record['owner'])


def _match_edge_runs(paths, neighbours):
    runs = [np.load(path, mmap_mode='r') for path in paths]
    try:
        previous = None
        owner = None
        paired = False
        for a, b, current in heapq.merge(*(_records(run) for run in runs)):
            edge = (a, b)
            if edge != previous:
                previous, owner, paired = edge, current, False
            elif paired:
                raise MeshTopologyError(f'Non-manifold edge {edge}: more than two triangles')
            else:
                neighbours.flat[owner] = current // 3
                neighbours.flat[current] = owner // 3
                paired = True
    finally:
        for run in runs:
            run._mmap.close()


def build_neighbours(elements, nvertex, vertices=None, chunk_size=DEFAULT_CHUNK_SIZE):
    """Return an int64 (N, 3) table, rejecting invalid triangle topology.

    chunk_size bounds the number of triangles sorted in memory at once. The
    triangle and resulting neighbour arrays still occupy memory; larger edge
    tables are written to temporary sorted runs under the system temporary dir.
    """
    triangles = integer_rows(elements, 'elem')
    if not len(triangles):
        raise MeshTopologyError('Mesh has no triangles')
    if not isinstance(chunk_size, (int, np.integer)) or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    if triangles.min() < 0 or triangles.max() >= nvertex:
        raise MeshTopologyError(f'elem contains a vertex index outside [0, {nvertex})')
    positions = None
    if vertices is not None:
        positions = np.asarray(vertices, dtype=float)
        if positions.ndim != 2 or positions.shape[0] != nvertex or positions.shape[1] not in (2, 3):
            raise MeshTopologyError('vertex must have shape (nvertex, 2) or (nvertex, 3)')
        if not np.isfinite(positions).all():
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

    neighbours = np.full(triangles.shape, -1, dtype=np.int64)
    if len(triangles) <= chunk_size:
        _match_sorted_edges(_edge_records(triangles, 0), neighbours)
    else:
        with tempfile.TemporaryDirectory(prefix='mesher-edges-') as directory:
            paths = []
            for start in range(0, len(triangles), chunk_size):
                path = os.path.join(directory, f'{start}.npy')
                np.save(path, _edge_records(triangles[start:start + chunk_size], start))
                paths.append(path)
            _match_edge_runs(paths, neighbours)

    for start in range(0, len(neighbours), chunk_size):
        block = neighbours[start:start + chunk_size]
        duplicates = ((block[:, 0] >= 0) & ((block[:, 0] == block[:, 1]) | (block[:, 0] == block[:, 2]))) | ((block[:, 1] >= 0) & (block[:, 1] == block[:, 2]))
        if duplicates.any():
            raise MeshTopologyError(f'Duplicate triangles involving face {start + np.flatnonzero(duplicates)[0]}')
    return neighbours


def validate_neighbours(values, expected=None):
    """Validate a neighbour graph and optionally its exact edge-slot mapping."""
    if len(values) == 0:
        raise MeshTopologyError('Neighbour table is empty. Run meshaudit --rebuild-neighbours on the mesh before partitioning')
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


def _boundary_defects(vertices, triangles, neighbours, tolerance):
    """Find coincident boundary vertices and T-junctions; do not weld them."""
    from scipy.spatial import cKDTree

    faces, slots = np.nonzero(neighbours == -1)
    a = triangles[faces, (slots + 1) % 3]
    b = triangles[faces, (slots + 2) % 3]
    ids = np.unique(np.concatenate((a, b)))
    if len(ids) < 2:
        return {'coincident_boundary_vertices': 0, 'hanging_nodes': 0, 'hanging_node_examples': []}
    xy = vertices[ids, :2]
    tree = cKDTree(xy)
    distances, _ = tree.query(xy, k=2)
    coincident = int(np.count_nonzero(distances[:, 1] <= tolerance))
    hanging = 0
    examples = []
    for face, va, vb in zip(faces, a, b):
        p, q = vertices[va, :2], vertices[vb, :2]
        delta = q - p
        length = np.linalg.norm(delta)
        candidates = tree.query_ball_point((p + q) / 2, length / 2 + tolerance)
        for candidate in candidates:
            vid = int(ids[candidate])
            if vid in (va, vb):
                continue
            offset = xy[candidate] - p
            along = np.dot(offset, delta) / length
            distance = abs(delta[0] * offset[1] - delta[1] * offset[0]) / length
            if tolerance < along < length - tolerance and distance <= tolerance:
                hanging += 1
                if len(examples) < 10:
                    examples.append({'vertex': vid, 'edge': [int(va), int(vb)], 'face': int(face)})
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
    report = {'errors': [], 'neighbour_errors': [], 'warnings': [],
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
    return report, neighbours
