import unittest

import numpy as np

from mesher.mesh_seams import reconcile_boundary_subdivisions
from mesher.mesh_topology import audit_mesh, build_neighbours, MeshTopologyError


def strip_mesh(left_y=(0, 1, 2), right_y=(0, .25, .5, .75, 1, 2)):
    """Two tiles with a partially matching seam: connected despite T-junctions."""
    vertices, triangles, ids = [], [], {}

    def vertex(x, y):
        key = (x, y)
        if key not in ids:
            ids[key] = len(vertices)
            vertices.append([x, y, 10.0])
        return ids[key]

    for x, ys in ((0, left_y), (1, right_y)):
        for lo, hi in zip(ys, ys[1:]):
            a, b, c, d = [vertex(*p) for p in ((x, lo), (x + 1, lo), (x + 1, hi), (x, hi))]
            triangles.extend([[a, b, c], [a, c, d]])
    return np.asarray(vertices), np.asarray(triangles)


def report(vertices, triangles, neighbours=None):
    if neighbours is None:
        neighbours = build_neighbours(triangles, len(vertices), vertices)
    return audit_mesh({'mesh': {'vertex': vertices.tolist(), 'elem': triangles.tolist(),
                               'neigh': neighbours.tolist(), 'nvertex': len(vertices),
                               'nelem': len(triangles)}})[0]


def areas(vertices, triangles):
    pts = vertices[triangles, :2]
    u, v = pts[:, 1] - pts[:, 0], pts[:, 2] - pts[:, 0]
    return .5 * (u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])


class MeshSeamTests(unittest.TestCase):
    def test_connected_mesh_three_nodes_on_one_edge(self):
        vertices, triangles = strip_mesh()
        original_vertices, original_triangles = vertices.copy(), triangles.copy()
        before = report(vertices, triangles)
        self.assertEqual(before['components'], 1)
        self.assertEqual(before['hanging_nodes'], 3)
        self.assertEqual(len({e['face'] for e in before['hanging_node_examples']}), 1)
        for chunk in (1, 1_000_000):
            result, neighbours, stats = reconcile_boundary_subdivisions(vertices, triangles, chunk_size=chunk)
            self.assertEqual(stats['split_faces'], 1)
            self.assertEqual(stats['added_triangles'], 3)
            after = report(vertices, result, neighbours)
            self.assertTrue(after['valid'], after)
            self.assertEqual(after['components'], 1)
            self.assertEqual(after['hanging_nodes'], 0)
            self.assertAlmostEqual(areas(vertices, result).sum(), areas(vertices, triangles).sum())
            self.assertTrue((areas(vertices, result) > 0).all())
        np.testing.assert_array_equal(vertices, original_vertices)
        np.testing.assert_array_equal(triangles, original_triangles)

    def test_interleaved_subdivisions_on_both_sides(self):
        vertices, triangles = strip_mesh((0, .4, 1), (0, .6, 1))
        for reverse in (False, True):
            source = triangles[:, ::-1] if reverse else triangles
            result, neighbours, stats = reconcile_boundary_subdivisions(vertices, source)
            self.assertEqual(stats['added_triangles'], 2)
            self.assertTrue(report(vertices, result, neighbours)['valid'])
            self.assertAlmostEqual(areas(vertices, result).sum(), areas(vertices, source).sum())
            self.assertTrue((areas(vertices, result) < 0).all() if reverse else (areas(vertices, result) > 0).all())

    def test_splits_on_all_three_edges_of_one_face(self):
        vertices = np.array([[0, 0, 0], [2, 0, 0], [0, 2, 0], [.5, 0, 0], [1.5, 0, 0],
                             [1, 1, 0], [0, 1, 0], [1, -2, 0], [2, 2, 0], [-2, 1, 0]])
        triangles = np.array([[0, 1, 2], [7, 0, 3], [7, 3, 4], [7, 4, 1],
                              [8, 1, 5], [8, 5, 2], [9, 2, 6], [9, 6, 0]])
        result, neighbours, stats = reconcile_boundary_subdivisions(vertices, triangles)
        self.assertEqual(stats['split_faces'], 1)
        self.assertEqual(stats['added_triangles'], 4)
        self.assertTrue(report(vertices, result, neighbours)['valid'])
        self.assertAlmostEqual(abs(areas(vertices, triangles)).sum(), abs(areas(vertices, result)).sum())

    def test_conforming_mesh_is_unchanged_and_repair_is_idempotent(self):
        vertices, triangles = strip_mesh()
        repaired, neighbours, _ = reconcile_boundary_subdivisions(vertices, triangles)
        again, again_neighbours, stats = reconcile_boundary_subdivisions(vertices, repaired)
        np.testing.assert_array_equal(again, repaired)
        np.testing.assert_array_equal(again_neighbours, neighbours)
        self.assertEqual(stats['added_triangles'], 0)

    def test_refuses_gap_and_overlap_without_mutating_input(self):
        vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, .5, 0], [2, .5, 0]])
        for triangles, message in (([[0, 1, 2], [4, 3, 1]], 'no matching'),
                                   ([[0, 1, 2], [4, 3, 1], [4, 2, 3]], 'overlapping')):
            points = vertices.copy()
            if message == 'overlapping':
                points[4] = [.75, .4, 0]
            elements = np.asarray(triangles)
            points_before, elements_before = points.copy(), elements.copy()
            with self.assertRaisesRegex(MeshTopologyError, message):
                reconcile_boundary_subdivisions(points, elements)
            np.testing.assert_array_equal(points, points_before)
            np.testing.assert_array_equal(elements, elements_before)

    def test_refuses_coincident_boundary_vertices(self):
        vertices, triangles = strip_mesh((0, 1), (0, 1))
        duplicate = len(vertices)
        vertices = np.vstack((vertices, vertices[1]))
        triangles[0, 1] = duplicate
        with self.assertRaisesRegex(MeshTopologyError, 'Coincident'):
            reconcile_boundary_subdivisions(vertices, triangles)

    def test_large_projected_coordinates_and_roundoff(self):
        vertices, triangles = strip_mesh()
        vertices[:, :2] += [500_000, 7_000_000]
        # Slight serialization-level deviation, within the default XY tolerance.
        vertices[(vertices[:, 0] == 500_001) & (vertices[:, 1] == 7_000_000.5), 0] += 1e-9
        result, neighbours, _ = reconcile_boundary_subdivisions(vertices, triangles)
        self.assertTrue(report(vertices, result, neighbours)['valid'])
        self.assertAlmostEqual(areas(vertices, result).sum(), areas(vertices, triangles).sum(), places=8)

    def test_repair_accepts_opposite_faces_with_different_winding(self):
        vertices, triangles = strip_mesh()
        triangles[4:] = triangles[4:, ::-1]
        result, neighbours, _ = reconcile_boundary_subdivisions(vertices, triangles)
        self.assertTrue(report(vertices, result, neighbours)['valid'])
        self.assertAlmostEqual(abs(areas(vertices, result)).sum(), abs(areas(vertices, triangles)).sum())

    def test_tolerance_must_be_valid(self):
        vertices, triangles = strip_mesh()
        for tolerance in (-1, float('nan'), float('inf')):
            with self.subTest(tolerance=tolerance), self.assertRaisesRegex(ValueError, 'tolerance'):
                reconcile_boundary_subdivisions(vertices, triangles, tolerance=tolerance)


if __name__ == '__main__':
    unittest.main()
