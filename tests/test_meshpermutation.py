import contextlib
import ctypes.util
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from mesher.mesh_topology import build_neighbours, MeshTopologyError
from mesher.tools import meshpermutation as permutation
from test_mesh_topology import square_mesh


class MeshPermutationTests(unittest.TestCase):
    def setUp(self):
        self.neighbours = build_neighbours(square_mesh()['mesh']['elem'], 4).tolist()
        self.plot_patch = patch.object(permutation, 'plot_mat_connectivity')
        self.plot_patch.start()
        self.addCleanup(self.plot_patch.stop)

    def test_invalid_graph_and_rank_count_never_call_metis(self):
        cases = [([], 2), (self.neighbours, 0), (self.neighbours, 1),
                 (self.neighbours, 3), ([[1, -1, -1], [-1, -1, -1]], 2)]
        with patch.object(permutation, '_load_metis') as loader:
            for rows, ranks in cases:
                with self.subTest(rows=rows, ranks=ranks):
                    with self.assertRaises(MeshTopologyError):
                        permutation.compute_metis_permutation_N(rows, ranks)
            loader.assert_not_called()

    def test_bad_metis_results_have_actionable_errors(self):
        neighbours = [[1, -1, -1], [0, 2, -1], [1, -1, -1]]
        for labels in ([], [0], [0, 1, 2], [-1, 0, 1], [0.0, 1.0, 0.0], [0, 0, 0]):
            with self.subTest(labels=labels):
                backend = SimpleNamespace(part_graph=Mock(return_value=(0, labels)))
                with patch.object(permutation, '_load_metis', return_value=backend):
                    with self.assertRaisesRegex(MeshTopologyError, 'METIS returned'):
                        permutation.compute_metis_permutation_N(neighbours, 2)

    def test_permutation_is_new_position_to_old_face(self):
        neighbours = [[1, -1, -1], [0, 2, -1], [1, -1, -1]]
        backend = SimpleNamespace(part_graph=Mock(return_value=(1, [1, 0, 1])))
        with patch.object(permutation, '_load_metis', return_value=backend), contextlib.redirect_stdout(io.StringIO()):
            sizes, order = permutation.compute_metis_permutation_N(neighbours, 2)
        self.assertEqual(sizes, [1, 2])
        self.assertEqual(order, [1, 0, 2])
        backend.part_graph.assert_called_once_with([[1], [0, 2], [1]], 2)

    def test_one_face_per_rank_has_an_exact_partition(self):
        with patch.object(permutation, '_load_metis') as loader, contextlib.redirect_stdout(io.StringIO()):
            sizes, order = permutation.compute_metis_permutation_N(self.neighbours, 2)
            loader.assert_not_called()
        self.assertEqual(sizes, [1, 1])
        self.assertEqual(order, [0, 1])

    def test_file_validation_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.mesh'
            output = Path(directory) / 'output.mesh'
            source.write_text(json.dumps(square_mesh()))
            output.write_text('untouched')
            args = {'infile': str(source), 'outfile': str(output), 'type': 'metis', 'nrank': 2}
            with patch.object(permutation, '_load_metis') as loader:
                with self.assertRaisesRegex(MeshTopologyError, 'Neighbour table is empty'):
                    permutation.append_global_cell_id_to_mesh_file(args)
                loader.assert_not_called()
            self.assertEqual(output.read_text(), 'untouched')

    def test_isolated_triangle_bandwidth_does_not_crash(self):
        matrix = permutation.convert_neighbor_list_to_compressed_matrix([[-1, -1, -1]], permutation.sparse.csr_matrix)
        self.assertEqual(permutation.compute_bandwidth_CSR(matrix), 0)

    def test_cli_rejects_empty_graph_before_loading_an_invalid_dll(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.mesh'
            source.write_text(json.dumps(square_mesh()))
            original = source.read_bytes()
            env = os.environ.copy()
            env['METIS_DLL'] = '/does-not-exist/libmetis.so'
            env['PYTHONPATH'] = str(Path(permutation.__file__).parents[2])
            result = subprocess.run([sys.executable, '-m', 'mesher.tools.meshpermutation',
                                     '-i', str(source), '-t', 'metis', '-n', '2'],
                                    env=env, cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn('Neighbour table is empty', result.stderr)
            self.assertNotIn('Traceback', result.stderr)
            self.assertEqual(source.read_bytes(), original)

    @unittest.skipUnless(importlib.util.find_spec('metis') and ctypes.util.find_library('metis'), 'METIS Python wrapper and shared library required')
    def test_real_metis_two_and_four_parts(self):
        width = 5
        triangles = []
        for y in range(width - 1):
            for x in range(width - 1):
                a = y * width + x
                triangles.extend([[a, a + 1, a + width + 1], [a, a + width + 1, a + width]])
        neighbours = build_neighbours(triangles, width * width).tolist()
        for ranks in (2, 4):
            with self.subTest(ranks=ranks), contextlib.redirect_stdout(io.StringIO()):
                sizes, order = permutation.compute_metis_permutation_N(neighbours, ranks)
                self.assertEqual(sum(sizes), len(triangles))
                self.assertTrue(all(n > 0 for n in sizes))
                self.assertEqual(sorted(order), list(range(len(triangles))))

    def test_repaired_mesh_file_partitions(self):
        from mesher.tools.meshaudit import main as audit_main

        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.mesh'
            repaired = Path(directory) / 'repaired.mesh'
            output = Path(directory) / 'partitioned.mesh'
            document = square_mesh()
            source.write_text(json.dumps(document))
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(audit_main([str(source), '--rebuild-neighbours', '-o', str(repaired)]), 0)
                permutation.append_global_cell_id_to_mesh_file({'infile': str(repaired), 'outfile': str(output), 'type': 'metis', 'nrank': 2})
            result = json.loads(output.read_text())['mesh']
            self.assertEqual(result['local_size'], [1, 1])
            self.assertEqual(sorted(result['cell_global_id']), [0, 1])
            self.assertEqual(result['elem'], document['mesh']['elem'])
            self.assertEqual(result['vertex'], document['mesh']['vertex'])
            self.assertEqual(result['version'], '2.0.0')

    @unittest.skipUnless(importlib.util.find_spec('metis') and ctypes.util.find_library('metis'), 'METIS Python wrapper and shared library required')
    def test_real_metis_after_seam_reconciliation(self):
        from mesher.mesh_seams import reconcile_boundary_subdivisions
        from test_mesh_seams import strip_mesh

        vertices, triangles = strip_mesh()
        repaired, neighbours, stats = reconcile_boundary_subdivisions(vertices, triangles)
        self.assertEqual(stats['added_triangles'], 3)
        for ranks in (2, 4):
            with self.subTest(ranks=ranks), contextlib.redirect_stdout(io.StringIO()):
                sizes, order = permutation.compute_metis_permutation_N(neighbours.tolist(), ranks)
                self.assertEqual(sum(sizes), len(repaired))
                self.assertTrue(all(n > 0 for n in sizes))
                self.assertEqual(sorted(order), list(range(len(repaired))))
