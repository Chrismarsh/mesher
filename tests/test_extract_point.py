import unittest

import numpy as np
from osgeo import gdal

from mesher.mesher_cli import build_mesh_from_arrays, extract_point


class ExtractPointTests(unittest.TestCase):
    def test_dem_vertices_and_nodata_fallback(self):
        nodata = -9999.0
        raster = gdal.GetDriverByName('MEM').Create('', 3, 3, 1, gdal.GDT_Float32)
        raster.SetGeoTransform((0, 1, 0, 3, 0, -1))
        band = raster.GetRasterBand(1)
        band.SetNoDataValue(nodata)
        band.WriteArray(np.array([[1, 2, 3], [4, nodata, 6], [7, 8, 9]], dtype='f4'))

        # One pixel west still has neighbours. Farther west the entire window is outside.
        self.assertEqual(extract_point(raster, -1.5, 1.5), 4.0)
        self.assertEqual(extract_point(raster, -3.5, 1.5), nodata)
        self.assertEqual(extract_point(raster, 1.5, 1.5), 5.0)  # NoData inside the DEM
        self.assertEqual(extract_point(raster, 1.5, 2.5), 2.0)
        self.assertEqual(extract_point(raster, 3, 0), 9.0)  # bottom-right boundary

        # MPI mesh assembly can then repair an unsampled vertex from its triangle.
        vertices = [(-3.5, 1.5), (1.5, 2.5), (2.5, 2.5)]
        mesh = build_mesh_from_arrays(vertices, [(0, 1, 2)], raster, band, None, True, False)
        self.assertEqual([v[2] for v in mesh['mesh']['vertex']], [2.5, 2.0, 3.0])
