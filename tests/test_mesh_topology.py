import unittest

import numpy as np

from mesher.mesh_topology import build_neighbours, MeshTopologyError, validate_neighbours


def square_mesh():
    return {'mesh': {'nvertex': 4, 'nelem': 2,
                     'vertex': [[0, 0, 10], [1, 0, 10], [1, 1, 10], [0, 1, 10]],
                     'elem': [[0, 1, 2], [0, 2, 3]], 'neigh': []}}


class MeshTopologyTests(unittest.TestCase):
    def test_opposite_vertex_slots_and_chunked_equivalence(self):
        mesh = square_mesh()['mesh']
        for chunk in (1, 10):
            with self.subTest(chunk=chunk):
                neighbours = build_neighbours(mesh['elem'], 4, mesh['vertex'], chunk)
                self.assertEqual(neighbours.tolist(), [[-1, 1, -1], [-1, -1, 0]])
                validate_neighbours(neighbours)

    def test_matching_tile_join(self):
        # Two adjacent square tiles after merging their shared vertex IDs.
        vertices = [[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1]]
        triangles = [[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]]
        neighbours = build_neighbours(triangles, 6, vertices, chunk_size=2)
        self.assertEqual(neighbours[0, 0], 3)
        self.assertEqual(neighbours[3, 1], 0)
        validate_neighbours(neighbours)

    def test_malformed_triangles_fail_in_both_sorting_paths(self):
        cases = [([[0, 1, 2], [2, 1, 0]], 3, 'Duplicate'),
                 ([[0, 1, 2], [1, 0, 3], [0, 1, 4]], 5, 'Non-manifold'),
                 ([[0, 0, 1]], 2, 'repeats'),
                 ([[0, 1, -1]], 3, 'outside'),
                 ([[0, 1, 3]], 3, 'outside'),
                 ([[0.0, 1.0, 2.0]], 3, 'integer'),
                 ([[0, 1]], 3, 'shape')]
        for elements, nvertex, message in cases:
            for chunk in (1, 100):
                with self.subTest(elements=elements, chunk=chunk):
                    with self.assertRaisesRegex(MeshTopologyError, message):
                        build_neighbours(elements, nvertex, chunk_size=chunk)

    def test_empty_and_degenerate_geometry(self):
        with self.assertRaisesRegex(MeshTopologyError, 'no triangles'):
            build_neighbours(np.empty((0, 3), dtype=int), 0)
        with self.assertRaisesRegex(MeshTopologyError, 'degenerate'):
            build_neighbours([[0, 1, 2]], 3, [[0, 0], [1, 1], [2, 2]])
        with self.assertRaisesRegex(MeshTopologyError, 'non-finite'):
            build_neighbours([[0, 1, 2]], 3, [[0, 0], [1, 1], [0, np.nan]])

    def test_graph_rejects_empty_out_of_range_self_and_asymmetric_rows(self):
        for rows, message in [([], 'empty'), ([[1, -1, -1]], 'indices'),
                              ([[0, -1, -1]], 'self-neighbour'),
                              ([[1, -1, -1], [-1, -1, -1]], 'reciprocal'),
                              ([[1, 1, -1], [0, -1, -1]], 'repeated')]:
            with self.subTest(rows=rows):
                with self.assertRaisesRegex(MeshTopologyError, message):
                    validate_neighbours(rows)

    def test_reciprocal_graph_must_also_have_correct_slots(self):
        expected = build_neighbours(square_mesh()['mesh']['elem'], 4)
        with self.assertRaisesRegex(MeshTopologyError, 'slots'):
            validate_neighbours([[1, -1, -1], [0, -1, -1]], expected)
        with self.assertRaisesRegex(MeshTopologyError, 'Expected 2'):
            validate_neighbours([[-1, -1, -1]], expected)

    def test_many_chunks_agree_with_single_sort(self):
        triangles = []
        width = 12
        for y in range(width - 1):
            for x in range(width - 1):
                a = y * width + x
                triangles.extend([[a, a + 1, a + width + 1], [a, a + width + 1, a + width]])
        expected = build_neighbours(triangles, width * width)
        result = build_neighbours(triangles, width * width, chunk_size=7)
        np.testing.assert_array_equal(result, expected)
        validate_neighbours(result)
