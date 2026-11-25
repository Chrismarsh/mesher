"""Restart through the CLI and both native meshing workers, including an interrupted restart."""
from contextlib import ExitStack, redirect_stdout
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from osgeo import gdal, osr

from mesher import mesher_cli as cli
from mesher.mesher_utls import mpi_resume
from mesher.mesher_utls.MPI_mesh_tiles import mesh_tile
from mesher.mesher_utls.MPI_merge_tiles import merge_tiles


class Walltime(Exception):
    pass


@unittest.skipUnless(os.environ.get('MESHER_TEST_EXE'), 'set MESHER_TEST_EXE to a built backend')
class RestartTests(unittest.TestCase):
    # Put the DEM, config and outputs in an isolated working directory so successive runs see
    # exactly the checkpoints written by that case.
    def setUp(self):
        old_cwd = os.getcwd()
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.addCleanup(os.chdir, old_cwd)
        self.root = Path(directory.name)
        self.config = self.root / 'mesh_config.py'
        self.output = self.root / 'mesh_config' / 'dem'
        self.plan = self.output / 'dem_mpi_checkpoint.json'
        self.dem = self.root / 'dem.tif'
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(32608)
        ds = gdal.GetDriverByName('GTiff').Create(str(self.dem), 48, 32, 1, gdal.GDT_Float32)
        ds.SetGeoTransform((500000, 1, 0, 7000000, 0, -1))
        ds.SetProjection(srs.ExportToWkt())
        ds.GetRasterBand(1).SetNoDataValue(-9999)
        ds.GetRasterBand(1).Fill(100)
        ds = None
        self.write_config()
        self.ran = []

    # Vary workers, tile count and tolerance independently to distinguish resource changes from
    # changes to mesh identity.
    def write_config(self, ranks=2, tiles=6, tolerance=1):
        self.config.write_text(f"""
dem_filename = 'dem.tif'
mesher_path = {os.environ['MESHER_TEST_EXE']!r}
mpi_mesh = True
MPI_nworkers = {ranks}
mpi_num_tiles = {tiles}
max_area = 32
min_area = 1
max_tolerance = {tolerance}
write_vtu = write_shp = False
""")
        # CLI invocations normally run in separate interpreters.
        sys.modules.pop('config', None)
        Path(importlib.util.cache_from_source(str(self.config))).unlink(missing_ok=True)

    # Exercise the normal CLI and checkpoint logic while controlling exactly which tile or pair
    # finishes before an interruption.
    def run_cli(self, *flags, stop=None):
        self.ran = []
        # Dispatch real worker functions locally so the test can stop between publications and
        # record which tasks a restart reruns.
        def run_tasks(tasks, launcher, workers, worker, label):
            for task in tasks:
                if worker == 'MPI_mesh_tiles.py':
                    mesh_tile(task)
                    name = task['tile_prefix']
                else:
                    merge_tiles(task)
                    name = Path(task['out_prefix']).name
                self.ran.append(name)
                if name == stop:
                    raise Walltime()

        with ExitStack() as stack:
            stack.enter_context(patch.object(sys, 'argv', ['mesher', *flags, str(self.config)]))
            stack.enter_context(patch.object(cli, '_run_mpi_mesh_tasks', run_tasks))
            # Parameter aggregation has its own MPI integration tests. The real
            # export still checks connectivity and samples elevations here.
            stack.enter_context(patch.object(cli, 'run_parameterization', return_value={}))
            if flags:
                stack.enter_context(patch.object(cli, 'project_dem', side_effect=AssertionError('reran warp')))
                stack.enter_context(patch.object(cli, 'prepare_domain_files', side_effect=AssertionError('reran polygonize')))
            log = stack.enter_context(redirect_stdout(io.StringIO()))
            cli.main()
            return log.getvalue()

    # Interrupt a merge stage, then replace a leaf. Only unfinished work and descendants of the
    # new leaf should need another solve.
    def test_restarts_mid_stage_and_after_replaced_leaf(self):
        first = 'dem_mpi_stage0_r0_c0'
        with self.assertRaises(Walltime):
            self.run_cli(stop=first)
        self.assertEqual(len(self.ran), 7)
        self.write_config(ranks=1, tiles=None)
        log = self.run_cli('--resume')
        self.assertIn('reusing 6/6, scheduling 0', log)
        self.assertIn('reusing 1/2, scheduling 1', log)
        self.assertEqual(len(self.ran), 4)
        layout = json.loads((self.output / 'dem_mpi_layout.json').read_text())
        final = self.output / layout['final_npz']
        with np.load(final) as arrays:
            expected = {k: arrays[k].copy() for k in arrays}
        self.assertTrue((self.output.parent / 'dem_mpi.mesh').is_file())

        # Simulate a truncated leaf, then a second interruption after replacing
        # it but before its old (still present) ancestors have been replaced.
        leaf = 'dem_mpi_tile_0'
        (self.output / (leaf + '.npz')).write_bytes(b'broken zip')
        with self.assertRaises(Walltime):
            self.run_cli('--resume', stop=leaf)
        self.assertEqual(self.ran, [leaf])
        self.assertFalse((self.output / 'dem_mpi_layout.json').exists())
        self.run_cli('--resume')
        self.assertEqual(self.ran, [first, 'dem_mpi_stage1_r0_c0', final.stem])
        self.run_cli('--resume')
        self.assertEqual(self.ran, [])
        with np.load(final) as recovered:
            for key in expected:
                np.testing.assert_array_equal(expected[key], recovered[key])
        self.write_config()
        self.run_cli()
        with np.load(final) as clean:
            for key in expected:
                np.testing.assert_array_equal(expected[key], clean[key])

    # Compare the evidence available in old tile markers with the stronger identity checks
    # possible once a run manifest exists.
    def test_legacy_adoption_checks_layout_and_tracked_resume_checks_inputs(self):
        first = 'dem_mpi_stage0_r0_c0'
        with self.assertRaises(Walltime):
            self.run_cli(stop=first)
        self.plan.unlink()
        for path in self.output.glob('dem_mpi_*.json'):
            marker = json.loads(path.read_text())
            marker.pop('checkpoint')
            path.write_text(json.dumps(marker))
        with self.assertRaisesRegex(ValueError, 'resume-legacy'):
            self.run_cli('--resume')
        self.write_config(tiles=4)
        with self.assertRaisesRegex(ValueError, 'legacy tile grid'):
            self.run_cli('--resume-legacy')
        self.assertFalse(self.plan.exists())
        self.write_config()
        self.run_cli('--resume-legacy')
        self.assertEqual(len(self.ran), 4)
        self.assertNotIn(first, self.ran)
        self.assertFalse(json.loads(self.plan.read_text())['legacy'])
        adopted = json.loads((self.output / (first + '.json')).read_text())['checkpoint']
        self.assertTrue(adopted['generation'].startswith('legacy:'))
        self.run_cli('--resume')
        self.assertEqual(self.ran, [])

        self.write_config(tolerance=30)
        with self.assertRaisesRegex(ValueError, 'criteria'):
            self.run_cli('--resume')
        self.assertEqual(self.ran, [])
        self.write_config()
        stat = self.dem.stat()
        os.utime(self.dem, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1000000))
        with self.assertRaisesRegex(ValueError, 'source inputs'):
            self.run_cli('--resume')
        self.assertEqual(self.ran, [])


class CheckpointTests(unittest.TestCase):
    # Interrupt the interval between archive and marker renames. A leftover NPZ cannot prove
    # that this task completed.
    def test_kill_between_publications_leaves_task_pending(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prefix = root/'one'
            staged = root/'staging.npz'
            task = dict(checkpoint_run='run', base_dir=directory + os.sep, tile_prefix='one',
                        tile_bbox=[0, 0, 10, 10], band_width=1)
            marker = dict(core_bbox=task['tile_bbox'], band_width=1, empty_tile=True, merge_method='remesh')
            plan = dict(run_id='run', legacy=False)
            # Recreate the worker's temporary archive before each attempted publication because
            # a successful rename consumes it.
            def stage():
                np.savez(staged, verts=np.empty((0, 3)), tris=np.empty((0, 3), dtype=int))
            stage()
            mpi_resume.publish(staged, prefix, marker.copy(), task)
            self.assertEqual(mpi_resume.unfinished([task], plan, 'tiles'), [])
            stage()
            replace = os.replace
            # Let the archive reach its final path and fail the marker rename to reproduce the
            # vulnerable publication interval.
            def interrupt_before_metadata(source, destination):
                if destination.suffix == '.json':
                    raise Walltime()
                replace(source, destination)
            with patch.object(os, 'replace', interrupt_before_metadata), self.assertRaises(Walltime):
                mpi_resume.publish(staged, prefix, marker.copy(), task)
            self.assertTrue(prefix.with_suffix('.npz').is_file())
            self.assertFalse(prefix.with_suffix('.json').exists())
            self.assertEqual(mpi_resume.unfinished([task], plan, 'tiles'), [task])

    # Reuse recorded warp paths and a config callback on resume so parameter restoration cannot
    # depend on relaunching GDAL workers.
    def test_projected_parameter_paths_and_classifier_are_restored_without_mpi(self):
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(32608)
        data = dict(soil=dict(file=['a.tif', 'b.tif'], method=['mode', 'mean'],
                              classifier=max, weight=0.25))
        prepared = dict(soil=['stored-a.tif', 'stored-b.tif'])
        with patch.object(cli, '_gdal_launch', side_effect=AssertionError('launched GDAL ranks')):
            total, weighted = cli.regularize_inputs('/cache', '', data, -1, 1, srs, 1, 0, 1, 0,
                False, None, 2, 'resume', prepared=prepared)
        self.assertEqual((total, weighted), (0.25, True))
        self.assertEqual(data['soil']['filename'], prepared['soil'])
        self.assertEqual(data['soil']['method'], ['mode', 'mean'])
        self.assertIsNone(data['soil']['file'])
        self.assertIsInstance(data['soil']['classifier'], bytes)

    # Competing writers to one output tree could mix checkpoint generations. Verify the second
    # lock acquisition fails.
    def test_concurrent_writer_is_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'dem'
            with mpi_resume.output_lock(output):
                with self.assertRaisesRegex(RuntimeError, 'Another Mesher'):
                    with mpi_resume.output_lock(output):
                        self.fail('acquired lock held by another writer')
            with mpi_resume.output_lock(output):
                pass
