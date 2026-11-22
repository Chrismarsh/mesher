import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from osgeo import gdal, osr

from mesher import mesher_cli as cli
from mesher.mesher_cli import _gdal_launch
from mesher.mesher_utls.MPI_regularize_inputs import regularize_inputs


class ProjectedParameterTests(unittest.TestCase):
    def setUp(self):
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.cache = Path(work.name) / 'projected_parameters'
        self.cache.mkdir()
        self.source = Path(work.name) / 'input.tif'
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(32613)
        ds = gdal.GetDriverByName('GTiff').Create(str(self.source), 4, 4, 1, gdal.GDT_Float32)
        ds.SetProjection(srs.ExportToWkt())
        ds.SetGeoTransform((0, 1, 0, 4, 0, -1))
        ds.GetRasterBand(1).Fill(12)
        ds = None

    def project(self, reuse=True, method='mean', xmax=4, fill=False, source=None):
        data = {'file': str(source or self.source), 'method': method}
        args = (str(self.cache), data, '', 'soil', -1, 1, 'EPSG:32613', xmax, 0, 4, 0, fill, reuse)
        return Path(regularize_inputs(args)['filename'][0])

    def test_reuse_after_default_run_and_refresh(self):
        raster = self.project(reuse=False)
        with patch.object(subprocess, 'check_call', side_effect=AssertionError('unexpected warp')):
            self.assertEqual(raster, self.project())

        self.assertNotEqual(raster, self.project(method='mode'))
        self.assertNotEqual(raster, self.project(xmax=3))
        self.assertNotEqual(raster, self.project(fill=True))
        raster.unlink()
        self.assertEqual(raster, self.project())
        self.assertTrue(raster.is_file())

        # The cache must also notice an updated tile behind a VRT.
        vrt = self.source.with_suffix('.vrt')
        gdal.BuildVRT(str(vrt), [str(self.source)]).FlushCache()
        mosaic = self.source.with_name('mosaic.vrt')
        gdal.BuildVRT(str(mosaic), [str(vrt)]).FlushCache()
        before = self.project(source=mosaic)
        stamp = self.source.stat()
        os.utime(self.source, ns=(stamp.st_atime_ns, stamp.st_mtime_ns + 1000000000))
        self.assertNotEqual(before, self.project(source=mosaic))

    @unittest.skipUnless(shutil.which('mpirun'), 'Requires an MPI launcher')
    def test_separate_gdal_job_runs_real_warps(self):
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(32613)
        inputs = {}
        for i in range(4):
            source = self.source.with_name(f'parameter {i}.tif')
            shutil.copyfile(self.source, source)
            inputs[str(i)] = {'file': str(source), 'method': 'mean'}
        # Record the cache setting received by each gdalwarp child.
        probe = self.cache / 'gdalwarp'
        probe.write_text('#!/bin/sh\necho "$GDAL_CACHEMAX" >> "$CACHE_LOG"\n'
                         'exec ' + shlex.quote(shutil.which('gdalwarp')) + ' "$@"\n')
        probe.chmod(0o755)
        log = self.cache / 'cache.log'
        launcher = 'mpirun -n {nworkers} --map-by node ' + shlex.quote(sys.executable)
        with patch.dict(os.environ, GDAL_CACHEMAX='64', CACHE_LOG=str(log)):
            cli.regularize_inputs(str(self.cache / 'shared cache'), str(self.cache) + '/',
                                  inputs, -1, 1, srs, 4, 0, 4, 0, False,
                                  'mesh-only-launcher', 4000, False, 2, launcher)
        self.assertEqual(log.read_text().splitlines(), ['64'] * 4)
        for data in inputs.values():
            self.assertIsNone(data['file'])
            self.assertEqual(data['method'], ['mean'])
            ds = gdal.Open(data['filename'][0])
            self.assertEqual(ds.ReadAsArray()[0, 0], 12)
        self.assertFalse(list((self.cache / 'shared cache').glob('.harmonize-*')))

    def test_failed_fill_never_publishes_an_unfilled_raster(self):
        with patch.object(gdal, 'FillNodata', side_effect=RuntimeError('failed fill')):
            with self.assertRaisesRegex(RuntimeError, 'failed fill'):
                self.project(fill=True)
        self.assertEqual(list(self.cache.iterdir()), [])
        self.assertTrue(self.project(fill=True).is_file())


class GdalLauncherTests(unittest.TestCase):
    site = ('--bind-to core --hetero-nodes --runtime-options fwd-environment '
            '--prtemca plm_ssh_args "-o SendEnv=PBS_JOBID -o LogLevel=ERROR" python')
    mesh = 'mpirun -n 4000 --map-by slot ' + site

    def test_gdal_placement_and_job_wide_cap(self):
        template = 'mpirun -n {nworkers} --map-by node ' + self.site
        self.assertEqual(_gdal_launch(self.mesh, 4000, template, 5, 40),
                         ('mpirun -n 5 --map-by node ' + self.site, 5))
        self.assertEqual(_gdal_launch(self.mesh, 4000, template, 5, 2),
                         ('mpirun -n 2 --map-by node ' + self.site, 2))
        self.assertEqual(_gdal_launch(self.mesh, 4000, None, 5, 40),
                         ('mpirun -n 5 --map-by slot ' + self.site, 5))
        self.assertEqual(_gdal_launch(None, 4000, None, 5, 2), (None, 2))

    def test_legacy_wrapper_and_fixed_gdal_launcher(self):
        self.assertEqual(_gdal_launch('submit-job', 4, None, 4, 2), ('submit-job', 4))
        with self.assertRaisesRegex(ValueError, 'MPI_gdal_exec_str'):
            _gdal_launch('submit-job', 4000, None, 5, 40)
        self.assertEqual(_gdal_launch('submit-job', 4000, 'submit-job {nworkers}', 5, 40),
                         ('submit-job 5', 5))
        self.assertEqual(_gdal_launch(self.mesh, 4000, 'mpiexec -np 2 python', 5, 40),
                         ('mpiexec -np 2 python', 2))

    def test_config_override_reaches_harmonization_only(self):
        # Stop at the first worker stage to exercise the CLI's actual config wiring.
        with tempfile.TemporaryDirectory() as directory, \
                patch.dict(os.environ, {'MESHER_NWORKERS_GDAL': '5'}), \
                patch.object(sys, 'argv', ['mesher', str(Path(directory) / 'case.py')]):
            self.addCleanup(os.chdir, os.getcwd())
            Path(sys.argv[1]).write_text(
                'dem_filename="input.tif"\nmax_area=100\nMPI_nworkers=4000\n'
                f'MPI_exec_str={self.mesh!r}\nnworkers_gdal=8\n'
                'MPI_gdal_exec_str="mpirun -n {nworkers} --map-by node python"\n'
                'parameter_files={"p":{"file":"input.tif","method":"mean"}}\n')
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(32613)
            source = gdal.GetDriverByName('GTiff').Create(str(Path(directory) / 'input.tif'), 2, 2, 1)
            source.SetProjection(srs.ExportToWkt())
            source.SetGeoTransform((0, 1, 0, 2, 0, -1))
            source = None
            class ReachedHarmonization(Exception):
                pass
            with patch.object(cli, 'find_gdal_prefix', return_value=''), \
                    patch.object(cli, 'regularize_inputs', side_effect=ReachedHarmonization) as launch:
                with self.assertRaises(ReachedHarmonization):
                    cli.main()
            mesh_launcher, mesh_workers, _, limit, gdal_launcher = launch.call_args.args[-5:]
            self.assertEqual((mesh_launcher, mesh_workers, limit), (self.mesh, 4000, 5))
            self.assertEqual(_gdal_launch(mesh_launcher, mesh_workers, gdal_launcher, limit, 40),
                             ('mpirun -n 5 --map-by node python', 5))
            os.environ['MESHER_NWORKERS_GDAL'] = '0'
            with self.assertRaisesRegex(ValueError, 'nworkers_gdal must be a positive integer'):
                cli.read_config(sys.argv[1])

    def test_no_parameters_launches_no_workers(self):
        with patch.object(subprocess, 'check_call', side_effect=AssertionError('unexpected launch')):
            self.assertEqual(cli.regularize_inputs(None, '', {}, 1, 1, None, 1, 0, 1, 0, False,
                                                   'opaque-wrapper', 4000, False, 5), (0, False))
