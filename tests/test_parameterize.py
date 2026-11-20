import io
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest
from contextlib import redirect_stdout

import cloudpickle
import numpy as np
from osgeo import gdal, ogr, osr

from mesher.mesher_utls import MPI_do_parameterize as parameterize


class ParameterizeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        self.config = self.work / 'config.py'
        self.config.write_text('mesher_debug = False\n')
        self.srs = osr.SpatialReference()
        self.srs.ImportFromEPSG(32608)
        self.gt = (0, 1, 0, 2, 0, -1)
        # Cover every pixel, including the diagonal edge at the far corner.
        self.mesh = {'mesh': {'vertex': [[-1, -1, 0], [5, -1, 0], [-1, 5, 0]],
                              'elem': [[0, 1, 2]], 'nelem': 1, 'nvertex': 3}}

    def raster(self, values, nodata=None, dtype=gdal.GDT_Float32, gt=None, name=None):
        values = np.asarray(values)
        driver = gdal.GetDriverByName('GTiff' if name else 'MEM')
        ds = driver.Create(str(self.work / name) if name else '', values.shape[1], values.shape[0], 1, dtype)
        ds.SetProjection(self.srs.ExportToWkt())
        ds.SetGeoTransform(self.gt if gt is None else gt)
        band = ds.GetRasterBand(1)
        band.WriteArray(values)
        if nodata is not None:
            band.SetNoDataValue(nodata)
        ds.FlushCache()
        return ds

    def aggregate(self, rasters, methods='mode', classifier=None, mesh=None):
        data = {'file': rasters, 'method': methods}
        if classifier is not None:
            data['classifier'] = cloudpickle.dumps(classifier)
        # Deliberately unrelated DEM grid: parameters must use their own grids.
        return parameterize.do_parameterize((100, 50, 0, 200, 0, -50), False, mesh or self.mesh,
            {'landcover': data}, 'landcover', {}, 1, 1, self.srs.ExportToProj4(), 0, str(self.config))['landcover']

    def test_nodata_is_excluded_before_mode_and_classifier(self):
        for nodata, dtype in ((-9999, gdal.GDT_Float32), (1, gdal.GDT_Float32),
                              (255, gdal.GDT_Byte), (np.nan, gdal.GDT_Float32)):
            with self.subTest(nodata=nodata, dtype=dtype):
                ds = self.raster([[4, nodata], [nodata, nodata]], nodata, dtype)
                self.assertEqual(self.aggregate([ds]), 4)
                # A typical catch-all class used to receive the NoData sentinel.
                self.assertEqual(self.aggregate([ds], classifier=lambda value: 2 if value == 4 else 18), 2)

    def test_zero_mask_scale_offset_and_missing_data(self):
        ds = self.raster([[0, 0], [4, 4]])
        self.assertEqual(self.aggregate([ds]), 0)  # Valid zero is not missing.
        band = ds.GetRasterBand(1)
        band.CreateMaskBand(0)
        band.GetMaskBand().WriteArray(np.array([[0, 0], [255, 255]], dtype=np.uint8))
        band.SetScale(.5)
        band.SetOffset(1)
        self.assertEqual(self.aggregate([ds], 'mean'), 3)
        band.SetScale(0)
        self.assertEqual(self.aggregate([ds], 'mean'), 1)
        empty = self.raster([[-9999, -9999], [-9999, -9999]], -9999)
        self.assertTrue(np.isnan(self.aggregate([empty])))
        # Keep the classifier's ability to explicitly handle a missing input.
        self.assertEqual(self.aggregate([empty], classifier=lambda value: 7 if np.isnan(value) else 0), 7)

    def test_multiple_raster_grids_and_callback_order(self):
        fine = self.raster([[2, 2], [2, -9999]], -9999)
        coarse = self.raster([[6]], gt=(0, 2, 0, 2, 0, -2))
        self.assertEqual(self.aggregate([fine, coarse], ['mean', 'max'],
                                       lambda a, b: 10*a+b), 26)
        self.assertEqual(self.aggregate([fine], lambda pixels: np.nanmax(pixels)+1,
                                       lambda value: 10*value), 30)

    def test_window_clips_all_sides_and_empty_intersections(self):
        offsets = parameterize.bbox_to_pixel_offsets
        self.assertEqual(offsets(self.gt, (-10, 1, 1, 10), 2, 2), (0, 0, 2, 2))
        self.assertEqual(offsets(self.gt, (1, 10, -10, 1), 2, 2), (1, 1, 1, 1))
        ds = self.raster([[4, 4], [4, 4]])
        for dx, dy in ((100, 0), (-100, 0), (0, 100), (0, -100)):
            with self.subTest(dx=dx, dy=dy):
                mesh = {'mesh': dict(self.mesh['mesh'], vertex=np.asarray(self.mesh['mesh']['vertex'], dtype=float)+[dx, dy, 0])}
                self.assertTrue(np.isnan(self.aggregate([ds], mesh=mesh)))

    def test_rotated_raster_uses_full_affine_transform(self):
        gt = (100, 1, .3, 200, .2, -1)
        ds = self.raster([[2, 2], [2, 2]], gt=gt)
        vertices = [[*gdal.ApplyGeoTransform(gt, x, y), 0] for x, y in ((-1, -1), (5, -1), (-1, 5))]
        mesh = {'mesh': dict(self.mesh['mesh'], vertex=vertices)}
        self.assertEqual(self.aggregate([ds], 'mean', mesh=mesh), 2)

    def run_worker(self, parameters, exact=True):
        args = [self.gt, False, self.mesh, parameters, {}, 2, 2, self.srs.ExportToProj4(), exact]
        path = self.work / 'args.pickle'
        path.write_bytes(cloudpickle.dumps(args))
        previous = os.getcwd()
        log = io.StringIO()
        try:
            os.chdir(self.work)
            with redirect_stdout(log):
                parameterize.main(str(path), False, str(self.config))
            with open('pickled_param_args_rets_0.pickle', 'rb') as stream:
                result = cloudpickle.load(stream)
        finally:
            os.chdir(previous)
        return result, log.getvalue()

    def test_exact_stats_and_legacy_callbacks_use_requested_methods(self):
        ds = self.raster([[1, 1], [4, 6]], name='values.tif')
        ds = None
        parameters = {method: {'filename': [str(self.work/'values.tif')], 'method': [method]}
                      for method in ('mean', 'average', 'mode', 'min', 'max')}
        parameters['custom'] = {'filename': [str(self.work/'values.tif')],
                                'method': [lambda pixels: np.nanmax(pixels)+100]}
        parameters['LAI'] = {'filename': [str(self.work/'values.tif')], 'method': ['mean'],
                             'classifier': cloudpickle.dumps(lambda value: value*2)}
        parameters['area'] = {}
        for exact in (False, True):
            with self.subTest(exact=exact):
                result, log = self.run_worker(parameters, exact)
                self.assertEqual(log, '')
                for key, value in {'mean': 3, 'average': 3, 'mode': 1, 'min': 1,
                                   'max': 6, 'custom': 106, 'LAI': 6}.items():
                    self.assertEqual(result[0][key], value, key)

    def test_unknown_exact_method_does_not_silently_compute_mean(self):
        ds = self.raster([[1, 1], [4, 6]], name='values.tif'); ds = None
        with self.assertRaisesRegex(ValueError, 'Unknown data aggregation method'):
            self.run_worker({'value': {'filename': [str(self.work/'values.tif')], 'method': ['typo']}})

    def test_configuration_is_loaded_once_for_all_triangles(self):
        self.config.write_text('from pathlib import Path\np = Path(__file__).with_suffix(".count")\n'
                               'p.write_text(str(int(p.read_text())+1) if p.exists() else "1")\n')
        ds = self.raster([[4, 4], [4, 4]])
        self.aggregate([ds])
        self.aggregate([ds])
        self.assertEqual(self.config.with_suffix('.count').read_text(), '1')

    @unittest.skipUnless(shutil.which('mpirun'), 'OpenMPI is required')
    def test_mpi_callbacks_from_config_and_empty_rank(self):
        self.config.write_text('import numpy as np\nmesher_debug = False\n'
                               'def classify(value):\n    return 2 if value == 4 else 18\n'
                               'def reduce(pixels):\n    return np.nanmax(pixels)+10\n')
        config = parameterize._load_config(str(self.config))
        ds = self.raster([[4, 4], [4, -9999]], -9999, name='classes.tif'); ds = None
        np.save(self.work/'vertices.npy', np.array([[0., 0, 0], [2., 0, 0], [0., 2, 0], [2., 2, 0]]))
        np.save(self.work/'triangles.npy', np.array([[0, 1, 2], [1, 3, 2]]))
        prep = {'gt': self.gt, 'is_geographic': False, 'RasterXSize': 2, 'RasterYSize': 2,
                'srs_proj4': self.srs.ExportToProj4(), 'use_exactextract': True,
                'verts_path': str(self.work/'vertices.npy'), 'elems_path': str(self.work/'triangles.npy'),
                'initial_conditions': {}, 'parameter_files': {
                    'landcover': {'filename': [str(self.work/'classes.tif')], 'method': ['mode'],
                                  'classifier': cloudpickle.dumps(config.classify)},
                    'raw': {'filename': [str(self.work/'classes.tif')], 'method': ['mode']},
                    'custom': {'filename': [str(self.work/'classes.tif')], 'method': [config.reduce]},
                    'area': {}}}
        inputs = self.work / 'prepare.pickle'
        inputs.write_bytes(cloudpickle.dumps(prep))
        command = ['mpirun', '--oversubscribe', '-n', '3', sys.executable, parameterize.__file__]
        for arguments in ([str(inputs), 'False', str(self.config), '--prepare-only'],
                          ['pickled_param_args_RANK.pickle', 'False', str(self.config)]):
            run = subprocess.run(command+arguments, cwd=self.work, text=True, capture_output=True, timeout=45)
            self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
            self.assertNotIn('per-pixel', run.stderr)
            self.assertEqual(run.stdout, '')
        results = [cloudpickle.loads((self.work/f'pickled_param_args_rets_{rank}.pickle').read_bytes())
                   for rank in range(3)]
        self.assertEqual(results[2], [])
        for result in results[:2]:
            self.assertEqual(result[0]['landcover'], 2)
            self.assertEqual(result[0]['raw'], 4)
            self.assertEqual(result[0]['custom'], 14)
