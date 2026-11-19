"""Reconcile subdivided triangle edges before MPI mesh parameter extraction.

This does not fill holes, weld vertices, or remesh arbitrary overlaps. A split
is accepted only if every resulting seam segment has an opposite-side face.
Inputs are never modified, including when validation fails.
"""

import numpy as np

from .mesh_topology import (MeshTopologyError, boundary_subdivisions,
                            build_neighbours, DEFAULT_CHUNK_SIZE, integer_rows)


def _cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def _area2(vertices, triangle):
    a, b, c = vertices[triangle, :2]
    return _cross(b - a, c - a)


def _split_face(triangle, splits):
    """Split original edges in turn, retaining the original winding."""
    children = [list(triangle)]
    for slot, nodes in sorted(splits.items()):
        a, b = triangle[(slot + 1) % 3], triangle[(slot + 2) % 3]
        for i, child in enumerate(children):
            if a in child and b in child:
                apex = next(v for v in child if v not in (a, b))
                chain = [a] + nodes + [b]
                children[i:i + 1] = [[apex, u, v] for u, v in zip(chain, chain[1:])]
                break
        else:
            raise MeshTopologyError(f'Cannot locate edge ({a}, {b}) while splitting a seam')
    return children


def reconcile_boundary_subdivisions(vertices, elements, tolerance=1e-8,
                                   chunk_size=DEFAULT_CHUNK_SIZE):
    """Return (triangles, neighbours, statistics) with conforming seam edges.

    Only existing vertices are used; coordinates and elevations are unchanged.
    This can change face count/order and must precede sampling face parameters.
    The tolerance is in mesh coordinate units. Local signed area is conserved
    within a tolerance-scaled bound; children must retain their parent's winding.
    This is a connectivity repair, not a guarantee of Delaunay/angle quality.
    """
    vertices = np.asarray(vertices, dtype=float)
    triangles = integer_rows(elements, 'elem')
    neighbours = build_neighbours(triangles, len(vertices), vertices, chunk_size)
    coincident, subdivisions = boundary_subdivisions(vertices, triangles, neighbours, tolerance)
    if coincident:
        raise MeshTopologyError('Coincident boundary vertices must be resolved before seam reconciliation')
    splits = {}
    required_edges = set()
    detections = 0
    for face, slot, nodes in subdivisions:
        splits.setdefault(face, {})[slot] = nodes
        detections += len(nodes)
        chain = [int(triangles[face, (slot + 1) % 3])] + nodes + [int(triangles[face, (slot + 2) % 3])]
        required_edges.update(tuple(sorted((a, b))) for a, b in zip(chain, chain[1:]))
    stats = {'split_faces': len(splits), 'hanging_node_detections': detections,
             'added_triangles': 0}
    if not splits:
        return triangles, neighbours, stats

    # Replace the original face with its first child, appending the others.
    # Keep only affected faces in Python objects; the full mesh stays in arrays.
    result = np.empty((len(triangles) + detections, 3), dtype=np.int64)
    result[:len(triangles)] = triangles
    cursor = len(triangles)
    changed_faces = []
    for face, edge_splits in sorted(splits.items()):
        children = _split_face(triangles[face], edge_splits)
        area = _area2(vertices, triangles[face])
        areas = np.array([_area2(vertices, child) for child in children])
        if np.any(np.sign(areas) != np.sign(area)):
            raise MeshTopologyError(f'Seam split in face {face} would invert or collapse a triangle')
        points = vertices[triangles[face], :2]
        perimeter = sum(np.linalg.norm(points[(i + 1) % 3] - points[i]) for i in range(3))
        scale = max(float(np.max(np.abs(points))), perimeter)
        roundoff = 64 * np.finfo(float).eps * scale * perimeter * len(children)
        if abs(float(areas.sum()) - area) > 2 * tolerance * perimeter + roundoff:
            raise MeshTopologyError(f'Seam split in face {face} would change its area beyond tolerance')
        result[face] = children[0]
        end = cursor + len(children) - 1
        result[cursor:end] = children[1:]
        changed_faces.extend([face] + list(range(cursor, end)))
        cursor = end
    if cursor != len(result):
        raise MeshTopologyError('Unexpected triangle count during seam reconciliation')

    # Building the graph also checks indices, duplicate/degenerate triangles and
    # non-manifold edges. Every split boundary segment must now be an interior
    # edge, with its two incident triangles lying on opposite sides in XY.
    neighbours = build_neighbours(result, len(vertices), vertices, chunk_size)
    checked = set()
    for face in changed_faces:
        tri = result[face]
        for slot in range(3):
            a, b = sorted((int(tri[(slot + 1) % 3]), int(tri[(slot + 2) % 3])))
            edge = (a, b)
            if edge not in required_edges or edge in checked:
                continue
            other = int(neighbours[face, slot])
            if other < 0:
                raise MeshTopologyError(f'Seam segment {edge} at face {face} has no matching opposite face; refusing to fill a gap')
            other_apex = next(v for v in result[other] if v not in edge)
            delta = vertices[b, :2] - vertices[a, :2]
            left = _cross(delta, vertices[tri[slot], :2] - vertices[a, :2])
            right = _cross(delta, vertices[other_apex, :2] - vertices[a, :2])
            if np.sign(left) == np.sign(right) or left == 0 or right == 0:
                raise MeshTopologyError(f'Seam segment {edge} has overlapping or degenerate incident faces')
            checked.add(edge)
    if checked != required_edges:
        raise MeshTopologyError('Not all subdivided seam segments were validated')
    coincident, remaining = boundary_subdivisions(vertices, result, neighbours, tolerance)
    if coincident or next(remaining, None) is not None:
        raise MeshTopologyError('Unresolved boundary defects after seam reconciliation')
    stats['added_triangles'] = len(result) - len(triangles)
    return result, neighbours, stats
