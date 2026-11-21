import unittest
from unittest.mock import patch

import numpy as np
from osgeo import gdal

from mesher.mesher_cli import build_mesh_from_arrays, extract_point
from mesher.mesh_topology import MeshTopologyError, audit_mesh


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
        np.testing.assert_array_equal(mesh['mesh']['neigh'], [[-1, -1, -1]])

    def test_mpi_array_export_includes_reciprocal_neighbours(self):
        raster = gdal.GetDriverByName('MEM').Create('', 3, 3, 1, gdal.GDT_Float32)
        raster.SetGeoTransform((0, 1, 0, 3, 0, -1))
        band = raster.GetRasterBand(1)
        band.SetNoDataValue(-9999)
        band.Fill(10)
        vertices = [(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)]
        mesh = build_mesh_from_arrays(vertices, [(0, 1, 2), (0, 2, 3)], raster, band, None, True, False)
        np.testing.assert_array_equal(mesh['mesh']['neigh'], [[-1, 1, -1], [-1, -1, 0]])
        self.assertEqual(mesh['mesh']['nelem'], 2)
        np.testing.assert_array_equal(mesh['mesh']['elem'], [[0, 1, 2], [0, 2, 3]])
        np.testing.assert_array_equal(np.asarray(mesh['mesh']['vertex'])[:, :2], vertices)
        self.assertTrue(audit_mesh(mesh)[0]['valid'])

    def test_hanging_nodes_rejected_before_sampling_without_modifying_inputs(self):
        vertices = np.array([[0, 0], [1, 0], [1, 1], [1, .5], [2, .5]], dtype=float)
        triangles = np.array([[0, 1, 2], [4, 3, 1], [4, 2, 3]])
        self.assert_export_rejects(vertices, triangles, 'Hanging boundary nodes')

    def test_coincident_boundary_vertices_rejected_before_sampling(self):
        vertices = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], dtype=float)
        triangles = np.array([[0, 1, 2], [4, 2, 3]])
        self.assert_export_rejects(vertices, triangles, 'Coincident boundary vertices')

    def assert_export_rejects(self, vertices, triangles, message):
        original_v, original_t = vertices.copy(), triangles.copy()
        with patch('mesher.mesher_cli.extract_point') as sample:
            with self.assertRaisesRegex(MeshTopologyError, message):
                build_mesh_from_arrays(vertices, triangles, None, None, None, True, False)
            sample.assert_not_called()
        np.testing.assert_array_equal(vertices, original_v)
        np.testing.assert_array_equal(triangles, original_t)
