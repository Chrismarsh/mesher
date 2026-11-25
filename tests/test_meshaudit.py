import contextlib
import copy
import io
import json
from pathlib import Path
import tempfile
import unittest

from mesher.mesh_topology import audit_mesh, build_neighbours
from mesher.tools.meshaudit import main
from test_mesh_topology import square_mesh


class MeshAuditTests(unittest.TestCase):
    # Keep the original and repaired mesh at different paths so each test can inspect what the
    # command wrote and check whether it touched its input.
    def setUp(self):
        self.work = tempfile.TemporaryDirectory()
        self.addCleanup(self.work.cleanup)
        self.source = Path(self.work.name) / 'input.mesh'
        self.output = Path(self.work.name) / 'repaired.mesh'

    # Invoke the audit CLI and collect its status and report so tests cover the same interface
    # used on saved meshes.
    def run_tool(self, document, *extra):
        self.source.write_text(json.dumps(document))
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            status = main([str(self.source), '--json'] + list(extra))
        return status, json.loads(stream.getvalue())

    # Exercise the empty-neighbour-table regression and verify that reconstruction keeps the
    # source file and the face geometry intact.
    def test_missing_neighbours_report_and_repair_preserve_input(self):
        original = square_mesh()
        original['mesh'].update(version='1.0.0', proj4='test', cell_global_id=[1, 0], local_size=[1, 1])
        original['metadata'] = {'keep': ['custom', 1.25]}
        status, report = self.run_tool(original)
        self.assertEqual(status, 1)
        self.assertEqual(report['components'], 1)
        self.assertIn('empty', report['neighbour_errors'][0])
        self.assertEqual(report['errors'], [])
        status, report = self.run_tool(original, '--rebuild-neighbours', '-o', str(self.output), '--chunk-size', '1')
        self.assertEqual(status, 0)
        self.assertTrue(report['valid'])
        expected = copy.deepcopy(original)
        expected['mesh']['neigh'] = [[-1, 1, -1], [-1, -1, 0]]
        self.assertEqual(json.loads(self.output.read_text()), expected)
        self.assertEqual(json.loads(self.source.read_text()), original)

    # Moving a neighbour into the wrong slot is repairable from connectivity even when all
    # referenced face IDs exist.
    def test_misplaced_slots_can_be_repaired(self):
        document = square_mesh()
        document['mesh']['neigh'] = [[1, -1, -1], [0, -1, -1]]
        status, report = self.run_tool(document, '--rebuild-neighbours', '-o', str(self.output))
        self.assertEqual(status, 0)
        self.assertIn('slots', report['previous_neighbour_errors'][0])

    # Changing neighbour rows cannot make an unsplit triangle meet a subdivided edge. Refuse to
    # publish a purported repair.
    def test_hanging_node_blocks_neighbour_only_repair(self):
        # Left tile has one interface edge, right tile has two at vertex 6.
        vertices = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0], [1, .5, 0]]
        triangles = [[0, 1, 2], [0, 2, 3], [6, 1, 4], [6, 4, 5], [6, 5, 2]]
        document = {'mesh': {'vertex': vertices, 'elem': triangles, 'neigh': [], 'nvertex': 7, 'nelem': 5}}
        status, report = self.run_tool(document, '--rebuild-neighbours', '-o', str(self.output))
        self.assertEqual(status, 1)
        self.assertEqual(report['components'], 2)
        self.assertEqual(report['hanging_nodes'], 1)
        self.assertEqual(report['hanging_node_examples'][0]['vertex'], 6)
        self.assertEqual(report['hanging_node_examples'][0]['vertex_xy'], [1, .5])
        self.assertFalse(report['neighbour_rebuild_allowed'])
        self.assertIn('blocked', report['recovery'])
        self.assertNotIn('Run meshaudit --rebuild-neighbours', report['neighbour_errors'][0])
        self.assertFalse(self.output.exists())

    # Duplicated boundary vertex IDs can disconnect geometrically touching faces. Check the
    # audit identifies the seam.
    def test_unwelded_seam_is_reported(self):
        document = square_mesh()
        document['mesh']['vertex'].extend([[0, 0, 10], [1, 1, 10]])
        document['mesh']['nvertex'] = 6
        document['mesh']['elem'][1] = [4, 5, 3]
        report, _ = audit_mesh(document)
        self.assertEqual(report['coincident_boundary_vertices'], 4)
        self.assertTrue(report['errors'])

    # Separate domain islands can be valid. A component count alone must not invalidate the
    # mesh.
    def test_intentionally_disconnected_components_are_warnings(self):
        vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [3, 0, 0], [4, 0, 0], [3, 1, 0]]
        document = {'mesh': {'vertex': vertices, 'elem': [[0, 1, 2], [3, 4, 5]],
                             'neigh': [[-1] * 3] * 2, 'nvertex': 6, 'nelem': 2}}
        report, _ = audit_mesh(document)
        self.assertTrue(report['valid'])
        self.assertEqual(report['components'], 2)
        self.assertEqual(report['largest_component_sizes'], [1, 1])
        self.assertTrue(report['warnings'])
        self.assertFalse(report['domain_coverage_checked'])

    # Supply invalid geometry and verify the repair path cannot publish an apparently successful
    # output file.
    def test_invalid_mesh_is_not_written(self):
        for change in ('count', 'index', 'finite', 'duplicate', 'degenerate'):
            document = square_mesh()
            m = document['mesh']
            if change == 'count': m['nelem'] = 999
            if change == 'index': m['elem'][0][0] = 99
            if change == 'finite': m['vertex'][0][2] = float('nan')
            if change == 'duplicate': m['elem'][1] = m['elem'][0][:]
            if change == 'degenerate': m['vertex'][2] = [2, 0, 10]
            with self.subTest(change=change):
                status, report = self.run_tool(document, '--rebuild-neighbours', '-o', str(self.output))
                self.assertEqual(status, 1)
                self.assertTrue(report['errors'])
                self.assertFalse(self.output.exists())

    # A repair destination must be new so a diagnostic run cannot destroy an input or an earlier
    # result.
    def test_refuses_input_and_existing_output(self):
        self.source.write_text(json.dumps(square_mesh()))
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main([str(self.source), '--rebuild-neighbours', '-o', str(self.source)])
        self.output.write_text('existing output')
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            main([str(self.source), '--rebuild-neighbours', '-o', str(self.output)])
        self.assertEqual(self.output.read_text(), 'existing output')
        self.assertEqual(sorted(p.name for p in Path(self.work.name).iterdir()), ['input.mesh', 'repaired.mesh'])
