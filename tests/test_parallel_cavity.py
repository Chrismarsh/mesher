"""Regressions for a large released seam, inter-piece conformance and PBS restart."""
import json
import os
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

from mpi4py import MPI
import numpy as np
from shapely.geometry import MultiPolygon, box

from mesher.mesher_utls import parallel_cavity as parallel
from mesher.mesher_utls import mpi_remesh
import test_remesh_tolerances as terrain


@unittest.skipUnless(os.environ.get('MESHER_TEST_EXE'), 'Set MESHER_TEST_EXE to the native backend')
class ParallelCavityTests(unittest.TestCase):
    def setUp(self):
        self.fixture = terrain.RemeshToleranceTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.options = dict(self.fixture.options, cavity_target_faces=10, max_tolerance=30,
                            _cavity_cache=self.fixture.work/'pieces')
        self.team = parallel.CavityTeam(MPI.COMM_SELF, self.options)
        self.options['_cavity_team'] = self.team

    def solve(self):
        return parallel.solve(self.fixture.domain, self.fixture.streams, self.options, 0, 40, 'test seam')

    def test_connected_seam_with_crossing_stream_at_all_tolerances(self):
        for tolerance in terrain.TOLERANCES_M:
            with self.subTest(tolerance=tolerance):
                self.options.update(max_tolerance=tolerance, _cavity_cache=self.fixture.work/str(tolerance))
                result = self.solve()
                self.assertGreater(result.stats['parallel_pieces'], 1)
                self.fixture.check_mesh(*result, tolerance)

    def test_restart_reuses_the_completed_piece_and_recovers_a_truncated_one(self):
        # Simulate PBS ending before gather or the pair checkpoint. Only two
        # ranks finished their first piece, and one of those files is damaged.
        def timeout(jobs):
            parallel._solve_piece(jobs[0], self.options)
            bad = parallel._solve_piece(jobs[1], self.options)[0]
            bad.write_bytes(b'truncated')
            raise TimeoutError
        with patch.object(self.team, 'map', side_effect=timeout), self.assertRaises(TimeoutError):
            self.solve()
        result = self.solve()
        self.assertEqual(result.stats['cached_piece_reads'], 1)
        self.fixture.check_mesh(*result, 30)
        with patch.object(mpi_remesh, 'remesh_domain', side_effect=AssertionError('unexpected native replay')):
            replay = self.solve()
        np.testing.assert_array_equal(result[1], replay[1])
        self.assertEqual(replay.stats['piece_solves'], 0)

    def test_killed_backend_is_fatal_to_the_distributed_solve(self):
        with patch.object(mpi_remesh, 'remesh_domain', side_effect=mpi_remesh.BackendTerminated('killed')), \
                self.assertRaisesRegex(mpi_remesh.BackendTerminated, 'killed'):
            self.solve()
        self.assertFalse(list(self.options['_cavity_cache'].glob('*.npz')))

    def test_aggregate_vertex_budget_applies_before_mesh_publication(self):
        # Cached pieces may individually fit a smaller bound while their
        # combined triangulation exceeds the single-cavity limit.
        self.solve()
        self.options['remesh_max_vertices'] = 150
        with self.assertRaisesRegex(mpi_remesh.PatchError, 'Patch vertex budget exceeded'):
            self.solve()

    def test_failed_artificial_cut_is_removed_locally(self):
        native = mpi_remesh.remesh_domain
        x, _, xmax, _ = self.fixture.domain.bounds
        troublesome_cut = x + (xmax-x)/4
        def reject_cut(domain, *args, **kwargs):
            if domain.bounds[2] == troublesome_cut:
                raise mpi_remesh.PatchError('infeasible cut')
            return native(domain, *args, **kwargs)
        with patch.object(mpi_remesh, 'remesh_domain', side_effect=reject_cut):
            result = self.solve()
        self.assertEqual(result.stats['parallel_pieces'], 3)
        self.fixture.check_mesh(*result, 30)

    def test_disconnected_components_and_hole_preserve_coverage(self):
        x, y, xmax, ymax = self.fixture.domain.bounds
        cut = x + (xmax-x)/2
        left = box(x, y, cut-15, ymax).difference(box(x+100, y+100, x+220, y+220))
        domain = MultiPolygon([left, box(cut+15, y, xmax, ymax)])
        result = parallel.solve(domain, [], self.options, 0, 40, 'test seam')
        self.assertGreater(result.stats['parallel_pieces'], 1)
        mpi_remesh.validate_patch(*result, domain)

    @unittest.skipUnless(shutil.which('mpirun'), 'MPI launcher required')
    def test_one_pending_pair_uses_multiple_ranks_and_publishes(self):
        from mesher.mesher_cli import run_mpi_merge_stage
        from mesher.mesher_utls.mpi_resume import archive_info
        from test_mpi_remesh import grid
        # A wide band releases all local faces. Use the actual launcher, task
        # manifest, worker teams, pair checkpoint and commit-marker ordering.
        f = self.fixture
        x, y, xmax, ymax = f.domain.bounds
        middle = (y+ymax)/2
        from osgeo import ogr
        shp = f.work/'domain.geojson'
        ds = ogr.GetDriverByName('GeoJSON').CreateDataSource(str(shp))
        layer = ds.CreateLayer('domain', geom_type=ogr.wkbPolygon)
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(ogr.CreateGeometryFromWkb(f.domain.wkb))
        layer.CreateFeature(feature)
        ds = None
        task = dict(f.options, max_tolerance=30, cavity_workers=3, cavity_target_faces=80,
                    remesh_band_width=2000, out_prefix=str(f.work/'pair'), constraints=[],
                    outer_polygon_shp=str(shp), checkpoint_run='test', mesher_debug=False)
        for name, bounds in [('a', (x, y, xmax, middle)), ('b', (x, middle, xmax, ymax))]:
            v, t = grid(bounds[0], bounds[2], bounds[1], bounds[3], step=75)
            archive, marker = f.work/(name+'.npz'), f.work/(name+'.json')
            np.savez(archive, verts=v, tris=t)
            marker.write_text(json.dumps(dict(core_bbox=bounds, merge_method='remesh',
                                              checkpoint={'generation': name})))
            task['tile_' + name] = {'npz': str(archive), 'meta': str(marker)}
        env = dict(OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
        with patch.dict(os.environ, env):
            run_mpi_merge_stage([task], 'mpirun --oversubscribe -n {nworkers} ' + sys.executable, 3)
        marker = json.loads((f.work/'pair.json').read_text())
        self.assertEqual(marker['checkpoint']['parents'], ['a', 'b'])
        self.assertEqual(marker['checkpoint']['archive'], archive_info(f.work/'pair.npz')[0])
        attempts = marker['remesh_report']['attempt_stats']
        self.assertTrue(any(a.get('parallel_pieces', 0) > 1 and a.get('cavity_ranks') == 3 for a in attempts))
        with np.load(f.work/'pair.npz') as data:
            mpi_remesh.validate_patch(data['verts'], data['tris'], f.domain)
