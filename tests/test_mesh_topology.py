import unittest

import numpy as np

from mesher.mesh_topology import build_neighbours, MeshTopologyError, validate_neighbours


# Two triangles sharing a diagonal give a minimal graph with one reciprocal neighbour and four
# boundary edges.
def square_mesh():
    return {'mesh': {'nvertex': 4, 'nelem': 2,
                     'vertex': [[0, 0, 10], [1, 0, 10], [1, 1, 10], [0, 1, 10]],
                     'elem': [[0, 1, 2], [0, 2, 3]], 'neigh': []}}


class MeshTopologyTests(unittest.TestCase):
    # Place the shared edge across a sort-run boundary and require both paths to use CHM's
    # opposite-vertex slot convention.
    def test_opposite_vertex_slots_and_chunked_equivalence(self):
        mesh = square_mesh()['mesh']
        for chunk in (1, 10):
            with self.subTest(chunk=chunk):
                neighbours = build_neighbours(mesh['elem'], 4, mesh['vertex'], chunk)
                self.assertEqual(neighbours.tolist(), [[-1, 1, -1], [-1, -1, 0]])
                validate_neighbours(neighbours)

    # A joined tile edge should appear as an internal neighbour pair even when it was a boundary
    # in each separate input.
    def test_matching_tile_join(self):
        # Two adjacent square tiles after merging their shared vertex IDs.
        vertices = [[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1]]
        triangles = [[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]]
        neighbours = build_neighbours(triangles, 6, vertices, chunk_size=2)
        self.assertEqual(neighbours[0, 0], 3)
        self.assertEqual(neighbours[3, 1], 0)
        validate_neighbours(neighbours)

    # Force the same invalid connectivity through the in-memory and external sorts so available
    # memory cannot change what Mesher accepts.
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

    # Distinguish empty or zero-area input from a valid face before its connectivity is used.
    def test_empty_and_degenerate_geometry(self):
        with self.assertRaisesRegex(MeshTopologyError, 'no triangles'):
            build_neighbours(np.empty((0, 3), dtype=int), 0)
        with self.assertRaisesRegex(MeshTopologyError, 'degenerate'):
            build_neighbours([[0, 1, 2]], 3, [[0, 0], [1, 1], [2, 2]])
        with self.assertRaisesRegex(MeshTopologyError, 'non-finite'):
            build_neighbours([[0, 1, 2]], 3, [[0, 0], [1, 1], [0, np.nan]])

    # Exercise invalid adjacency encodings before they reach a partitioner that assumes a
    # well-formed graph.
    def test_graph_rejects_empty_out_of_range_self_and_asymmetric_rows(self):
        for rows, message in [([], 'empty'), ([[1, -1, -1]], 'indices'),
                              ([[0, -1, -1]], 'self-neighbour'),
                              ([[1, -1, -1], [-1, -1, -1]], 'reciprocal'),
                              ([[1, 1, -1], [0, -1, -1]], 'repeated')]:
            with self.subTest(rows=rows):
                with self.assertRaisesRegex(MeshTopologyError, message):
                    validate_neighbours(rows)

    # A symmetric graph can still put the neighbour across the wrong triangle edge. Compare the
    # slots to the actual vertices.
    def test_reciprocal_graph_must_also_have_correct_slots(self):
        expected = build_neighbours(square_mesh()['mesh']['elem'], 4)
        with self.assertRaisesRegex(MeshTopologyError, 'slots'):
            validate_neighbours([[1, -1, -1], [0, -1, -1]], expected)
        with self.assertRaisesRegex(MeshTopologyError, 'Expected 2'):
            validate_neighbours([[-1, -1, -1]], expected)

    # Force enough tiny runs to exercise the multi-pass merge and compare neighbour IDs to a
    # single global sort.
    def test_many_chunks_agree_with_single_sort(self):
        triangles = []
        width = 12
        for y in range(width - 1):
            for x in range(width - 1):
                a = y * width + x
                triangles.extend([[a, a + 1, a + width + 1], [a, a + width + 1, a + width]])
        expected = build_neighbours(triangles, width * width)
        # More than 64 runs exercises bounded fan-in across multiple passes.
        result = build_neighbours(triangles, width * width, chunk_size=2)
        np.testing.assert_array_equal(result, expected)
        validate_neighbours(result)

    def test_interleaved_runs_match_opposite_slots_with_64bit_vertex_ids(self):
        # Every triangle meets at the centre, but its two radial neighbours can
        # land in any sort run. Compute their row IDs directly from the fan.
        count = 385
        base = 2**40
        order = np.random.default_rng(27).permutation(count)
        row = np.argsort(order)
        triangles = np.column_stack((np.full(count, base), base+1+order, base+1+(order+1) % count))
        expected = np.column_stack((np.full(count, -1), row[(order+1) % count], row[(order-1) % count]))
        for chunk in (1, 5, np.int8(64)):
            with self.subTest(chunk=chunk):
                actual = build_neighbours(triangles, base+count+1, chunk_size=chunk)
                np.testing.assert_array_equal(actual, expected)

    def test_third_owner_in_a_distant_run_is_rejected(self):
        # Put the three owners into separate generations of the external merge.
        # A boundary pair must not be finalized before the third record is seen.
        triangles = np.arange(3, 3+3*193).reshape(-1, 3)
        triangles[[0, 80, 192], :2] = (0, 1)
        with self.assertRaisesRegex(MeshTopologyError, r'Non-manifold edge \(0, 1\)'):
            build_neighbours(triangles, int(triangles.max()+1), chunk_size=1)
