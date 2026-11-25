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
from collections import namedtuple

from mesher.mesh_topology import build_neighbours, MeshTopologyError
from mesher.tools import meshpermutation as permutation
from test_mesh_topology import square_mesh


# Control METIS's returned labels so tests can exercise validation and permutation direction
# independently of the installed library.
def metis_stub(labels):
    return SimpleNamespace(part_graph=Mock(return_value=(1, labels)), idx_t=ctypes.c_int32,
                           METIS_Graph=namedtuple('Graph', 'nvtxs ncon xadj adjncy vwgt vsize adjwgt'))


class MeshPermutationTests(unittest.TestCase):
    # Use two faces sharing a known edge and disable plot output so the ordering tests can
    # inspect adjacency without writing diagnostic images.
    def setUp(self):
        self.neighbours = build_neighbours(square_mesh()['mesh']['elem'], 4).tolist()
        self.plot_patch = patch.object(permutation, 'plot_mat_connectivity')
        self.plot_patch.start()
        self.addCleanup(self.plot_patch.stop)

    # An invalid rank request or graph should fail in Python before the native partitioner can
    # receive it.
    def test_invalid_graph_and_rank_count_never_call_metis(self):
        cases = [([], 2), (self.neighbours, 0), (self.neighbours, 1),
                 (self.neighbours, 3), ([[1, -1, -1], [-1, -1, -1]], 2)]
        with patch.object(permutation, '_load_metis') as loader:
            for rows, ranks in cases:
                with self.subTest(rows=rows, ranks=ranks):
                    with self.assertRaises(MeshTopologyError):
                        permutation.compute_metis_permutation_N(rows, ranks)
            loader.assert_not_called()

    # Return malformed labels and empty parts from the stub to exercise validation of native
    # partitioner output.
    def test_bad_metis_results_have_actionable_errors(self):
        neighbours = [[1, -1, -1], [0, 2, -1], [1, -1, -1]]
        for labels in ([], [0], [0, 1, 2], [-1, 0, 1], [0.0, 1.0, 0.0], [0, 0, 0]):
            with self.subTest(labels=labels):
                backend = metis_stub(labels)
                with patch.object(permutation, '_load_metis', return_value=backend):
                    with self.assertRaisesRegex(MeshTopologyError, 'METIS returned'):
                        permutation.compute_metis_permutation_N(neighbours, 2)

    # Use an ordering whose inverse differs so the direction of the mapping and the order within
    # each rank are observable.
    def test_permutation_is_new_position_to_old_face(self):
        neighbours = [[1, -1, -1], [0, 2, -1], [1, -1, -1]]
        backend = metis_stub([1, 0, 1])
        with patch.object(permutation, '_load_metis', return_value=backend), contextlib.redirect_stdout(io.StringIO()):
            sizes, order = permutation.compute_metis_permutation_N(neighbours, 2)
        self.assertEqual(sizes, [1, 2])
        self.assertEqual(order, [1, 0, 2])
        graph_arg, parts = backend.part_graph.call_args[0]
        self.assertEqual(parts, 2)
        self.assertEqual(graph_arg.xadj[:4], [0, 1, 3, 4])
        self.assertEqual(graph_arg.adjncy[:4], [1, 0, 2, 1])

    # With as many ranks as faces the assignment is determined. Some METIS versions otherwise
    # return an empty part.
    def test_one_face_per_rank_has_an_exact_partition(self):
        with patch.object(permutation, '_load_metis') as loader, contextlib.redirect_stdout(io.StringIO()):
            sizes, order = permutation.compute_metis_permutation_N(self.neighbours, 2)
            loader.assert_not_called()
        self.assertEqual(sizes, [1, 1])
        self.assertEqual(order, [0, 1])

    # A rejected input must leave the requested destination intact, even when that file already
    # exists.
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

    # A one-face adjacency matrix has no off-diagonal edges but still needs a defined bandwidth
    # diagnostic.
    def test_isolated_triangle_bandwidth_does_not_crash(self):
        matrix = permutation.convert_neighbor_list_to_compressed_matrix([[-1, -1, -1]], permutation.sparse.csr_matrix)
        self.assertEqual(permutation.compute_bandwidth_CSR(matrix), 0)

    # Select an unusable METIS library to verify invalid connectivity is diagnosed before native
    # library discovery.
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

    # Use the installed partitioner to check valid, nonempty rank blocks for more than one
    # partition count.
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

    # Pass a reconstructed neighbour table into the normal permutation workflow to check the
    # audit and partition tools agree.
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
