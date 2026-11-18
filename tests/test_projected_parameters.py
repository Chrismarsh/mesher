import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from osgeo import gdal, osr

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

    def test_failed_fill_never_publishes_an_unfilled_raster(self):
        with patch.object(gdal, 'FillNodata', side_effect=RuntimeError('failed fill')):
            with self.assertRaisesRegex(RuntimeError, 'failed fill'):
                self.project(fill=True)
        self.assertEqual(list(self.cache.iterdir()), [])
        self.assertTrue(self.project(fill=True).is_file())
