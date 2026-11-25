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
from unittest.mock import Mock, patch
from types import SimpleNamespace

import cloudpickle
import numpy as np
from osgeo import gdal, ogr, osr

from mesher.mesher_utls import MPI_do_parameterize as parameterize


class ParameterizeTests(unittest.TestCase):
    # Provide a small mesh, raster reference system and isolated config so reducers can be
    # exercised with controlled pixel values.
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

    # Build a raster with chosen samples and georeferencing to test masking, scale and alignment
    # independently of input preparation.
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

    # Evaluate one parameter for a single test face through the GDAL path, including its
    # optional classifier.
    def aggregate(self, rasters, methods='mode', classifier=None, mesh=None):
        data = {'file': rasters, 'method': methods}
        if classifier is not None:
            data['classifier'] = cloudpickle.dumps(classifier)
        # Deliberately unrelated DEM grid: parameters must use their own grids.
        return parameterize.do_parameterize((100, 50, 0, 200, 0, -50), False, mesh or self.mesh,
            {'landcover': data}, 'landcover', {}, 1, 1, self.srs.ExportToProj4(), 0, str(self.config))['landcover']

    # A frequent NoData sentinel must not win the categorical mode or appear as a genuine class
    # in the callback.
    def test_nodata_is_excluded_before_mode_and_classifier(self):
        for nodata, dtype in ((-9999, gdal.GDT_Float32), (1, gdal.GDT_Float32),
                              (255, gdal.GDT_Byte), (np.nan, gdal.GDT_Float32)):
            with self.subTest(nodata=nodata, dtype=dtype):
                ds = self.raster([[4, nodata], [nodata, nodata]], nodata, dtype)
                self.assertEqual(self.aggregate([ds]), 4)
                # A typical catch-all class used to receive the NoData sentinel.
                self.assertEqual(self.aggregate([ds], classifier=lambda value: 2 if value == 4 else 18), 2)

    # Distinguish valid zero from masked pixels and verify scale and offset are applied only
    # after raw NoData has been identified.
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

    # Use rasters on different grids so each window must be computed separately and their
    # aggregates still reach the classifier in config order.
    def test_multiple_raster_grids_and_callback_order(self):
        fine = self.raster([[2, 2], [2, -9999]], -9999)
        coarse = self.raster([[6]], gt=(0, 2, 0, 2, 0, -2))
        self.assertEqual(self.aggregate([fine, coarse], ['mean', 'max'],
                                       lambda a, b: 10*a+b), 26)
        self.assertEqual(self.aggregate([fine], lambda pixels: np.nanmax(pixels)+1,
                                       lambda value: 10*value), 30)

    # Exercise clips on each raster edge as well as an empty intersection to catch negative or
    # oversized GDAL windows.
    def test_window_clips_all_sides_and_empty_intersections(self):
        offsets = parameterize.bbox_to_pixel_offsets
        self.assertEqual(offsets(self.gt, (-10, 1, 1, 10), 2, 2), (0, 0, 2, 2))
        self.assertEqual(offsets(self.gt, (1, 10, -10, 1), 2, 2), (1, 1, 1, 1))
        ds = self.raster([[4, 4], [4, 4]])
        for dx, dy in ((100, 0), (-100, 0), (0, 100), (0, -100)):
            with self.subTest(dx=dx, dy=dy):
                mesh = {'mesh': dict(self.mesh['mesh'], vertex=np.asarray(self.mesh['mesh']['vertex'], dtype=float)+[dx, dy, 0])}
                self.assertTrue(np.isnan(self.aggregate([ds], mesh=mesh)))

    # A rotated input makes a north-up shortcut give the wrong sampling window even if it works
    # on the usual DEM grid.
    def test_rotated_raster_uses_full_affine_transform(self):
        gt = (100, 1, .3, 200, .2, -1)
        ds = self.raster([[2, 2], [2, 2]], gt=gt)
        vertices = [[*gdal.ApplyGeoTransform(gt, x, y), 0] for x, y in ((-1, -1), (5, -1), (-1, 5))]
        mesh = {'mesh': dict(self.mesh['mesh'], vertex=vertices)}
        self.assertEqual(self.aggregate([ds], 'mean', mesh=mesh), 2)

    # Build the manifest and invoke the parameter worker so tests see the same row ordering,
    # raster choices and callbacks used by the launcher.
    def run_worker(self, parameters, exact=True, geographic=False):
        np.save(self.work/'vertices.npy', np.asarray(self.mesh['mesh']['vertex'], dtype=float))
        np.save(self.work/'elements.npy', np.asarray(self.mesh['mesh']['elem'], dtype=int).reshape(-1, 3))
        args = {'gt': self.gt, 'is_geographic': geographic, 'parameter_files': parameters,
                'initial_conditions': {}, 'RasterXSize': 2, 'RasterYSize': 2,
                'srs_proj4': self.srs.ExportToProj4(), 'use_exactextract': exact,
                'verts_path': str(self.work/'vertices.npy'), 'elems_path': str(self.work/'elements.npy'),
                'nranks': 1, 'output_dir': str(self.work)}
        path = self.work / 'args.pickle'
        path.write_bytes(cloudpickle.dumps(args))
        previous = os.getcwd()
        log = io.StringIO()
        try:
            os.chdir(self.work)
            with redirect_stdout(log):
                parameterize.main(str(path), False, str(self.config))
            with open(self.work/'result_0.pickle', 'rb') as stream:
                result = cloudpickle.load(stream)['triangles']
        finally:
            os.chdir(previous)
        return result, log.getvalue()

    # Combine built-in statistics with a Python reducer to verify each parameter uses its
    # configured sampling path.
    def test_exact_stats_classifiers_and_custom_reducers_use_requested_methods(self):
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

    # Use partially covered cells whose area-weighted result differs from an equal-pixel
    # average, then inspect what the classifier received.
    def test_classifiers_use_coverage_weighted_aggregates(self):
        # Pixel overlap areas are 5/6, 1/2, 1/6. Class 9 wins by area,
        # although class 1 wins by touched-pixel count. The numeric mean
        # is 10/3 by area, versus 6 with the legacy sampler.
        self.mesh['mesh']['vertex'] = [[0., 0., 0.], [3., 0., 0.], [0., 1., 0.]]
        for name, values, gt in (
                ('classes.tif', [[9, 1, 1]], (0, 1, 0, 1, 0, -1)),
                ('numeric.tif', [[0, 6, 12]], (0, 1, 0, 1, 0, -1)),
                ('coarse.tif', [[7]], (0, 3, 0, 1, 0, -1))):
            ds = self.raster(values, gt=gt, name=name); ds = None
        parameters = {
            'watermask': {'filename': str(self.work/'classes.tif'), 'method': 'mode',
                          'classifier': cloudpickle.dumps(lambda value: value)},
            'LAI': {'filename': [str(self.work/'numeric.tif')], 'method': ['mean'],
                    # A nonlinear callback proves classification occurs after
                    # averaging, rather than averaging classified pixels.
                    'classifier': cloudpickle.dumps(lambda value: value**2)},
            'combined': {'filename': [str(self.work/'classes.tif'), str(self.work/'coarse.tif')],
                         'method': ['mode', 'mean'],
                         'classifier': cloudpickle.dumps(lambda water, canopy: 10*water+canopy)},
            'custom': {'filename': [str(self.work/'classes.tif'), str(self.work/'numeric.tif')],
                       'method': ['mode', lambda pixels: np.nanmean(pixels)],
                       'classifier': cloudpickle.dumps(lambda a, b: 10*a+b)}}
        result, log = self.run_worker(parameters)
        self.assertEqual(log, '')
        self.assertEqual(result[0]['watermask'], 9)
        self.assertAlmostEqual(result[0]['LAI'], (10/3)**2, places=5)
        self.assertEqual(result[0]['combined'], 97)
        # If any reducer is callable, all inputs keep legacy weighting.
        self.assertEqual(result[0]['custom'], 16)
        parameters['watermask']['filename'] = [str(self.work/'classes.tif')]
        result, _ = self.run_worker(parameters, exact=False)
        self.assertEqual(result[0]['watermask'], 1)
        self.assertEqual(result[0]['LAI'], 36)
        self.assertEqual(result[0]['combined'], 17)

    # Choose a sentinel whose scaled value looks ordinary so a missing cell cannot silently
    # enter the classifier.
    def test_exact_classifier_masks_raw_missing_values_before_scaling(self):
        for nodata, dtype in ((-9999, gdal.GDT_Float32), (1, gdal.GDT_Float32),
                              (255, gdal.GDT_Byte), (np.nan, gdal.GDT_Float32)):
            with self.subTest(nodata=nodata):
                ds = self.raster([[4, nodata], [nodata, nodata]], nodata, dtype, name='valid.tif')
                band = ds.GetRasterBand(1)
                # An explicit mask must not override the raw NoData exclusion.
                band.CreateMaskBand(0)
                band.GetMaskBand().WriteArray(np.full((2, 2), 255, dtype=np.uint8))
                band.SetScale(.5); band.SetOffset(1)
                band = None; ds = None
                parameters = {method: {'filename': [str(self.work/'valid.tif')], 'method': [method],
                                       'classifier': cloudpickle.dumps(lambda value: value*10)}
                              for method in ('mode', 'mean', 'average', 'min', 'max')}
                result, _ = self.run_worker(parameters)
                for method in parameters:
                    self.assertEqual(result[0][method], 30)

        ds = self.raster([[0, np.inf], [np.nan, -np.inf]], name='finite.tif'); ds = None
        parameters = {method: {'filename': [str(self.work/'finite.tif')], 'method': [method],
                               'classifier': cloudpickle.dumps(lambda value: value)}
                      for method in ('mode', 'mean', 'min', 'max')}
        result, _ = self.run_worker(parameters)
        for method in parameters:
            self.assertEqual(result[0][method], 0)

        for scale, expected in ((.5, 3), (0, 1)):
            ds = self.raster([[0, 0], [4, 4]], name='masked.tif')
            band = ds.GetRasterBand(1)
            band.CreateMaskBand(0)
            band.GetMaskBand().WriteArray(np.array([[0, 0], [255, 255]], dtype=np.uint8))
            band.SetScale(scale); band.SetOffset(1)
            band = None; ds = None
            result, _ = self.run_worker({'value': {'filename': [str(self.work/'masked.tif')],
                'method': ['mean'], 'classifier': cloudpickle.dumps(lambda value: value)}})
            self.assertEqual(result[0]['value'], expected)

    # Exercise empty raster intersections and classifier return conventions so missing face
    # values survive result collection.
    def test_exact_missing_aggregates_and_classifier_outputs(self):
        for name, values in (('missing.tif', [[-9999, -9999], [-9999, -9999]]),
                             ('valid.tif', [[4, 4], [4, 4]])):
            ds = self.raster(values, -9999, name=name); ds = None
        parameters = {}
        for method in ('mode', 'mean', 'average', 'min', 'max'):
            data = {'filename': [str(self.work/'missing.tif')], 'method': [method]}
            parameters[method] = data
            parameters[method+'_classified'] = dict(data,
                classifier=cloudpickle.dumps(lambda value: 7 if np.isnan(value) else 18))
        parameters['partial'] = {'filename': [str(self.work/'missing.tif'), str(self.work/'valid.tif')],
            'method': ['mode', 'mean'], 'classifier': cloudpickle.dumps(lambda a, b: b if np.isnan(a) else a)}
        result, _ = self.run_worker(parameters)
        for method in ('mode', 'mean', 'average', 'min', 'max'):
            self.assertTrue(np.isnan(result[0][method]))
            self.assertEqual(result[0][method+'_classified'], 7)
        self.assertEqual(result[0]['partial'], 4)

        # No overlap has the same missing-data contract as all-NoData pixels.
        self.mesh['mesh']['vertex'] = (np.array(self.mesh['mesh']['vertex'], dtype=float)+[100, 100, 0]).tolist()
        result, _ = self.run_worker(parameters)
        self.assertTrue(np.isnan(result[0]['partial']))
        for method in ('mode', 'mean', 'average', 'min', 'max'):
            self.assertTrue(np.isnan(result[0][method]))
            self.assertEqual(result[0][method+'_classified'], 7)
        data = {'filename': [str(self.work/'valid.tif')], 'method': ['mean']}
        for output in (np.nan, -9999):
            result, _ = self.run_worker({'value': dict(data, classifier=cloudpickle.dumps(lambda value: output))})
            self.assertTrue(np.isnan(result[0]['value']))
        with self.assertRaisesRegex(ValueError, 'classifier function for value returned None'):
            self.run_worker({'value': dict(data, classifier=cloudpickle.dumps(lambda value: None))})

    # Count calls on a multi-pixel face. Classification happens after zonal aggregation and the
    # diagnostic should identify that path.
    def test_exact_classifier_runs_once_per_triangle_and_logs_backend(self):
        self.config.write_text('mesher_debug = True\n')
        ds = self.raster([[4, 4], [4, 4]], name='values.tif'); ds = None
        self.mesh['mesh'].update(vertex=[[0., 0., 0.], [2., 0., 0.], [0., 2., 0.], [2., 2., 0.]],
                                 elem=[[0, 1, 2], [1, 3, 2]], nelem=2, nvertex=4)
        calls = self.work/'classifier_calls.txt'

        # Leave evidence of each invocation so the test can distinguish one aggregate per face
        # from a call per pixel.
        def classify(value):
            with calls.open('a') as stream:
                stream.write(f'{value}\n')
            return value+1

        result, log = self.run_worker({'value': {'filename': [str(self.work/'values.tif')],
            'method': ['mean'], 'classifier': cloudpickle.dumps(classify)}})
        self.assertEqual(calls.read_text().splitlines(), ['4.0', '4.0'])
        self.assertEqual([row['value'] for row in result], [5, 5])
        self.assertIn('exactextract aggregation', log)
        self.assertNotIn('legacy ALL_TOUCHED', log)

    # An unrecognized aggregation method must fail explicitly because substituting a mean can
    # produce plausible but wrong parameters.
    def test_unknown_exact_method_does_not_silently_compute_mean(self):
        ds = self.raster([[1, 1], [4, 6]], name='values.tif'); ds = None
        with self.assertRaisesRegex(ValueError, 'Unknown data aggregation method'):
            self.run_worker({'value': {'filename': [str(self.work/'values.tif')], 'method': ['typo']}})

    # Parameters sharing a source can reuse one scan, but a subsequent invocation must read the
    # source again.
    def test_shared_raster_statistics_are_batched_and_reused_only_within_run(self):
        for name, values in (('shared.tif', [[1, 1], [4, 6]]), ('other.tif', [[10, 10], [10, 10]])):
            ds = self.raster(values, name=name); ds = None
        shared, other = str(self.work/'shared.tif'), str(self.work/'other.tif')
        parameters = {
            'mean': {'filename': [shared], 'method': ['mean']},
            'other': {'filename': [other], 'method': ['mean']},
            'combined': {'filename': [shared, other, shared], 'method': ['average', 'max', 'mode'],
                         'classifier': cloudpickle.dumps(lambda a, b, c: 100*a+10*b+c)},
            'repeat': {'filename': [shared, shared], 'method': ['mean', 'average'],
                       'classifier': cloudpickle.dumps(lambda a, b: a+b)}}
        original = parameterize.exact_extract
        calls = []

        # Record which statistics each extraction requested so the test can distinguish a shared
        # scan from repeated per-parameter reads.
        def extract(source, features, stats):
            calls.append((source.ds.GetDescription(), stats))
            return original(source, features, stats)

        with patch.object(parameterize, 'exact_extract', extract):
            result, _ = self.run_worker(parameters)
            self.assertEqual(result[0]['mean'], 3)
            self.assertEqual(result[0]['combined'], 401)
            self.assertEqual(result[0]['repeat'], 6)
            self.assertEqual(calls, [(shared, ['mean', 'majority']), (other, ['mean', 'max'])])
            ds = self.raster([[2, 2], [2, 2]], name='shared.tif'); ds = None
            result, _ = self.run_worker(parameters)
            self.assertEqual(result[0]['mean'], 2)
            self.assertEqual(result[0]['combined'], 302)
            self.assertEqual(result[0]['repeat'], 4)
            self.assertEqual(len(calls), 4)

    # Compare the batched area calculation with OGR on large projected values where subtracting
    # global-coordinate products loses precision.
    def test_projected_areas_match_ogr_for_large_coordinates_and_both_orientations(self):
        # Include a narrow face and a degenerate face; chunk boundaries must
        # not affect face ordering or the sign of the unsigned area.
        vertices = np.array([[0., 0., 0.], [10., 0., 9.], [0., 8., -1.],
                             [10., .0001, 3.], [5., 0., 0.]])+[-1900000, 3500000, 0]
        triangles = np.array([[0, 1, 2], [2, 1, 0], [0, 1, 3], [0, 1, 4]])
        self.mesh['mesh'].update(vertex=vertices, elem=triangles, nelem=4, nvertex=5)
        expected = []
        for i in range(4):
            expected.append(parameterize.do_parameterize(self.gt, False, self.mesh, {'area': {}},
                'area', {}, 2, 2, self.srs.ExportToProj4(), i, str(self.config))['area'])
        areas = parameterize._projected_triangle_areas(self.mesh, chunk_size=2)
        np.testing.assert_allclose(areas, expected, rtol=1e-12)
        self.assertEqual(areas[0], 40)
        self.assertEqual(areas[1], 40)
        self.assertEqual(areas[3], 0)
        # Neither projected worker path needs OGR per-triangle area objects.
        with patch.object(parameterize, 'do_parameterize', side_effect=AssertionError('unexpected OGR call')):
            for exact in (False, True):
                result, _ = self.run_worker({'area': {}}, exact)
                np.testing.assert_allclose([row['area'] for row in result], expected, rtol=1e-12)
        self.mesh['mesh'].update(elem=np.empty((0, 3), dtype=int), nelem=0)
        result, _ = self.run_worker({'area': {}})
        self.assertEqual(result, [])

    # Geographic faces need area in square metres. Check the path using an equal-area projection
    # rather than planar degree units.
    def test_geographic_areas_keep_existing_projected_geometry_calculation(self):
        self.srs = osr.SpatialReference(); self.srs.ImportFromEPSG(4326)
        self.mesh['mesh']['vertex'] = [[-135., 68., 0.], [-134.99, 68., 0.], [-135., 68.01, 0.]]
        expected = parameterize.do_parameterize(self.gt, True, self.mesh, {'area': {}}, 'area', {},
            2, 2, self.srs.ExportToProj4(), 0, str(self.config))['area']
        self.assertGreater(expected, 0)
        with patch.object(parameterize, '_projected_triangle_areas', side_effect=AssertionError('unprojected area')):
            for exact in (False, True):
                result, _ = self.run_worker({'area': {}}, exact, geographic=True)
                self.assertEqual(result[0]['area'], expected)

    # A config can have side effects. Multiple triangles in one worker should share the same
    # loaded module.
    def test_configuration_is_loaded_once_for_all_triangles(self):
        self.config.write_text('from pathlib import Path\np = Path(__file__).with_suffix(".count")\n'
                               'p.write_text(str(int(p.read_text())+1) if p.exists() else "1")\n')
        ds = self.raster([[4, 4], [4, 4]])
        self.aggregate([ds])
        self.aggregate([ds])
        self.assertEqual(self.config.with_suffix('.count').read_text(), '1')

    # Check disjoint face ranges and the backing mappings so transferring work to more ranks
    # does not duplicate the prepared mesh.
    def test_rank_ownership_is_contiguous_and_arrays_remain_file_backed(self):
        np.save(self.work/'vertices.npy', np.zeros((100000, 3)))
        elements = np.arange(30).reshape(10, 3)
        np.save(self.work/'elements.npy', elements)
        prep = {'verts_path': str(self.work/'vertices.npy'), 'elems_path': str(self.work/'elements.npy')}
        # Only paths are serialized, regardless of the global vertex count.
        self.assertLess(len(cloudpickle.dumps(prep)), 1000)
        for nranks in (1, 4, 16):
            expected = np.array_split(elements, nranks)
            previous = 0
            for rank in range(nranks):
                mesh, start, stop = parameterize._load_rank_mesh(prep, rank, nranks)
                self.assertEqual(start, previous)
                self.assertEqual(stop-start, len(expected[rank]))
                previous = stop
                for key in ('vertex', 'elem'):
                    if mesh['mesh'][key].size:
                        self.assertIsInstance(mesh['mesh'][key], np.memmap)
                        self.assertIsNotNone(mesh['mesh'][key].filename)
                    else:
                        self.assertEqual(mesh['mesh'][key].nbytes, 0)
                    self.assertFalse(mesh['mesh'][key].flags.writeable)
                np.testing.assert_array_equal(mesh['mesh']['elem'], expected[rank])
            self.assertEqual(previous, len(elements))

    # The manifest's row ownership assumes a specific communicator size. Catch a launcher
    # mismatch before writing results.
    def test_worker_rejects_mismatched_launcher_size(self):
        inputs = self.work/'mismatched.pickle'
        inputs.write_bytes(cloudpickle.dumps({'nranks': 2}))
        with self.assertRaisesRegex(ValueError, 'expected 2 MPI ranks'):
            parameterize.main(str(inputs), False, str(self.config))

    # Simulate missing or inconsistent worker results so the parent cannot silently accept an
    # incomplete parameter column.
    def test_launcher_propagates_failure_and_rejects_invalid_rank_results(self):
        from mesher import mesher_cli
        metadata = {'gt': self.gt, 'is_geographic': False, 'parameter_files': {},
                    'initial_conditions': {}, 'RasterXSize': 2, 'RasterYSize': 2,
                    'srs_proj4': self.srs.ExportToProj4(), 'use_exactextract': True}
        # Exercise the Spawn path independently of a cluster launcher. Failed
        # jobs cannot reuse unrelated pickles in the working directory.
        for mode in ('failure', 'missing', 'wrong_ids'):
            # Stand in for MPI.Spawn to inspect the shared manifest and control the result files
            # seen by the parent.
            def spawn(*args, **kwargs):
                inputs = Path(kwargs['args'][1])
                prep = cloudpickle.loads(inputs.read_bytes())
                self.assertNotIn('mesh', prep)
                self.assertNotIn('vertex', prep)
                if mode == 'failure':
                    raise RuntimeError('worker failed')
                if mode == 'wrong_ids':
                    (inputs.parent/'result_0.pickle').write_bytes(cloudpickle.dumps(
                        {'rank': 0, 'start': 0, 'stop': 1, 'triangles': [{'id': 8}]}))
                return Mock()

            fake_mpi = SimpleNamespace(COMM_SELF=SimpleNamespace(Spawn=Mock(side_effect=spawn)))
            with patch.object(mesher_cli, 'MPI', fake_mpi), self.subTest(mode=mode):
                error = {'failure': RuntimeError, 'missing': FileNotFoundError, 'wrong_ids': ValueError}[mode]
                with self.assertRaises(error):
                    mesher_cli.run_parameterization(self.mesh, metadata, str(self.work), None, 1, str(self.config))
                self.assertEqual(fake_mpi.COMM_SELF.Spawn.call_count, 1)
                self.assertEqual(list(self.work.glob('mesher-parameters-*')), [])

    # Run config-defined callbacks across real Python processes, including a rank with no faces
    # that must still complete.
    @unittest.skipUnless(shutil.which('mpirun'), 'OpenMPI is required')
    def test_mpi_callbacks_from_config_and_empty_rank(self):
        self.config.write_text('import numpy as np\nmesher_debug = False\n'
                               'def classify(water, numeric):\n    return 100*water+numeric\n'
                               'def reduce(pixels):\n    return np.nanmax(pixels)+10\n')
        config = parameterize._load_config(str(self.config))
        for name, values in (('classes.tif', [[9, 1, 1]]), ('numeric.tif', [[0, 6, 12]])):
            ds = self.raster(values, gt=(0, 1, 0, 1, 0, -1), name=name); ds = None
        np.save(self.work/'vertices.npy', np.array([[0., 0, 0], [3., 0, 0], [0., 1, 0], [3., 1, 0]]))
        np.save(self.work/'triangles.npy', np.array([[0, 1, 2], [1, 3, 2]]))
        prep = {'gt': self.gt, 'is_geographic': False, 'RasterXSize': 2, 'RasterYSize': 2,
                'srs_proj4': self.srs.ExportToProj4(), 'use_exactextract': True,
                'verts_path': str(self.work/'vertices.npy'), 'elems_path': str(self.work/'triangles.npy'),
                'nranks': 3, 'output_dir': str(self.work),
                'initial_conditions': {}, 'parameter_files': {
                    'landcover': {'filename': [str(self.work/'classes.tif'), str(self.work/'numeric.tif')],
                                  'method': ['mode', 'mean'],
                                  'classifier': cloudpickle.dumps(config.classify)},
                    'raw': {'filename': [str(self.work/'classes.tif')], 'method': ['mode']},
                    'custom': {'filename': [str(self.work/'classes.tif')], 'method': [config.reduce]},
                    'area': {}}}
        inputs = self.work / 'prepare.pickle'
        inputs.write_bytes(cloudpickle.dumps(prep))
        command = ['mpirun', '--oversubscribe', '-n', '3', sys.executable, parameterize.__file__]
        run = subprocess.run(command+[str(inputs), 'False', str(self.config)],
                             cwd=self.work, text=True, capture_output=True, timeout=45)
        self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
        self.assertNotIn('per-pixel', run.stderr)
        self.assertEqual(run.stdout, '')
        batches = [cloudpickle.loads((self.work/f'result_{rank}.pickle').read_bytes()) for rank in range(3)]
        results = [batch['triangles'] for batch in batches]
        self.assertEqual(results[2], [])
        self.assertTrue(inputs.exists())  # The parent owns shared-input cleanup.
        self.assertEqual(list(self.work.glob('pickled_param_args_*.pickle')), [])
        for result, water, numeric in zip(results[:2], (9, 1), (10/3, 26/3)):
            self.assertAlmostEqual(result[0]['landcover'], 100*water+numeric, places=5)
            self.assertEqual(result[0]['raw'], water)
            self.assertEqual(result[0]['custom'], 19)

        # The CLI launcher uses the same single-pass protocol. Verify quoting,
        # global IDs, exact rank collection and temporary-file cleanup.
        from mesher import mesher_cli
        import shlex
        root = self.work/'shared directory with spaces'; root.mkdir()
        config_with_spaces = self.work/'config with spaces.py'
        config_with_spaces.write_text(self.config.read_text())
        mesh = {'mesh': {'vertex': np.load(prep['verts_path']), 'elem': np.load(prep['elems_path'])}}
        command = f'mpirun --oversubscribe -n 3 {shlex.quote(sys.executable)}'
        with redirect_stdout(io.StringIO()):
            collected = mesher_cli.run_parameterization(mesh, prep, str(root), command, 3, str(config_with_spaces))
        self.assertEqual(collected, results[0]+results[1])
        self.assertEqual([row['id'] for row in collected], [0, 1])
        self.assertEqual(list(root.iterdir()), [])
        with tempfile.TemporaryDirectory(dir=root) as output_dir:
            with redirect_stdout(io.StringIO()):
                columns = mesher_cli.run_parameterization(mesh, dict(prep, batch_size=1), str(root),
                    command.replace('-n 3', '-n {nworkers}'), 3, str(config_with_spaces), output_dir=output_dir)
            for key, column in columns.items():
                self.assertIsInstance(column, np.memmap)
                np.testing.assert_allclose(column, [row[key] for row in collected])
        self.assertEqual(list(root.iterdir()), [])

        spawn_script = '''
import cloudpickle, numpy as np, sys
from mesher import mesher_cli
from mesher.mesher_utls import MPI_do_parameterize as worker
worker._load_config(sys.argv[2])
with open(sys.argv[1], 'rb') as stream:
    prep = cloudpickle.load(stream)
mesh = {'mesh': {'vertex': np.load(prep['verts_path']), 'elem': np.load(prep['elems_path'])}}
result = mesher_cli.run_parameterization(mesh, prep, sys.argv[3], None, 3, sys.argv[2])
with open(sys.argv[4], 'wb') as stream:
    cloudpickle.dump(result, stream)
'''
        spawned = self.work/'spawned.pickle'
        run = subprocess.run([sys.executable, '-c', spawn_script, str(inputs), str(config_with_spaces),
                              str(root), str(spawned)], cwd=self.work, text=True,
                             capture_output=True, timeout=45)
        self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
        self.assertEqual(cloudpickle.loads(spawned.read_bytes()), collected)
        self.assertEqual(list(root.iterdir()), [])
