"""Exercise elevation tolerances through the native solver and both seam merge paths."""
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
from osgeo import gdal, ogr, osr
from shapely.geometry import LineString, Polygon, box

from mesher.mesher_utls.disk_merge import merge_disk_cavity
from mesher.mesher_utls.mpi_remesh import merge_cavity, remesh_domain, tile_domain, validate_patch
from test_mpi_remesh import audit


TOLERANCES_M = (0.25, 0.75, 1, 5, 20, 30, 100, 300, 500)


def triangle_errors(vertices, triangles, elevation, transform, srs, min_area):
    """Independently fit vertex planes in DEM pixel coordinates over ALL_TOUCHED cells."""
    x0, dx, _, y0, _, dy = transform
    vectors = ogr.GetDriverByName('Memory').CreateDataSource('')
    layer = vectors.CreateLayer('face', srs, ogr.wkbPolygon)
    feature = ogr.Feature(layer.GetLayerDefn())
    masks = gdal.GetDriverByName('MEM')
    errors = []
    for xy in vertices[triangles, :2]:
        polygon = Polygon(xy)
        if polygon.area <= min_area:
            continue  # The area floor exempts elevation error only.
        pixels = np.floor((xy - (x0, y0)) / (dx, dy)).astype(int)
        lo, hi = pixels.min(axis=0), pixels.max(axis=0)
        local = pixels - lo
        matrix = np.column_stack((local, np.ones(3)))
        heights = elevation[pixels[:, 1], pixels[:, 0]]
        plane = np.linalg.solve(matrix, heights)

        cols, rows = (hi - lo + 1).tolist()
        mask = masks.Create('', cols, rows, 1, gdal.GDT_Byte)
        mask.SetGeoTransform((x0 + lo[0]*dx, dx, 0, y0 + lo[1]*dy, 0, dy))
        mask.SetProjection(srs.ExportToWkt())
        feature.SetGeometry(ogr.CreateGeometryFromWkb(polygon.wkb))
        layer.CreateFeature(feature)
        gdal.RasterizeLayer(mask, [1], layer, burn_values=[1], options=['ALL_TOUCHED=TRUE'])
        layer.DeleteFeature(feature.GetFID())
        feature.SetFID(-1)
        j, i = np.nonzero(mask.ReadAsArray())
        observed = elevation[j + lo[1], i + lo[0]]
        residual = observed - (plane[0]*i + plane[1]*j + plane[2])
        errors.append(np.sqrt(np.mean(residual**2)))
    return np.asarray(errors)


@unittest.skipUnless(os.environ.get('MESHER_TEST_EXE'), 'Set MESHER_TEST_EXE to the rebuilt meshercxx')
class RemeshToleranceTests(unittest.TestCase):
    # Construct projected terrain with a known elevation field and a physical constraint so the
    # same features can be checked at every tolerance.
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.work = Path(temp.name)
        x, y = -1817361., 3677022.
        self.domain = box(x, y, x+1350, y+1050)
        self.streams = [LineString([(x, y+90), (x+1350, y+950)])]
        # Off-grid boundaries avoid ambiguity about which pixel owns a vertex.
        self.transform = (x-51., 8., 0., y+1103., 0., -8.)
        self.srs = osr.SpatialReference()
        self.srs.ImportFromEPSG(32608)
        row, col = np.indices((146, 184))
        # A mountain on the seam and lower relief away from it exercise
        # different parts of the error range on the same source DEM.
        self.elevation = (1700*np.exp(-((col-90)**2+(row-70)**2)/200.)
                          + 42*np.exp(-((col-35)**2+(row-105)**2)/210.)
                          - 13*np.exp(-((col-138)**2+(row-43)**2)/130.)).astype(np.float32)
        dem = self.work/'terrain.tif'
        ds = gdal.GetDriverByName('GTiff').Create(str(dem), 184, 146, 1, gdal.GDT_Float32)
        ds.SetGeoTransform(self.transform)
        ds.SetProjection(self.srs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(self.elevation)
        ds = None
        self.options = dict(mesher_path=os.environ['MESHER_TEST_EXE'], dem_path=str(dem),
                            max_area=60000, min_area=1200, errormetric='rmse',
                            remesh_band_width=150, remesh_max_attempts=6, remesh_max_vertices=10000)

    # Independently measure terrain error and check topology and stream coverage so the solver
    # cannot validate itself using its own quality report.
    def check_mesh(self, vertices, triangles, tolerance):
        validate_patch(vertices, triangles, self.domain, self.streams)
        topology = audit(vertices, triangles)
        self.assertTrue(topology['valid'], topology)
        self.assertEqual(topology['components'], 1)
        xy = vertices[triangles, :2]
        for corner in range(3):
            a = xy[:, (corner+1) % 3] - xy[:, corner]
            b = xy[:, (corner+2) % 3] - xy[:, corner]
            twice_area = np.abs(a[:, 0]*b[:, 1] - a[:, 1]*b[:, 0])
            self.assertGreaterEqual(np.arctan2(twice_area, (a*b).sum(axis=1)).min(),
                                    np.arcsin(np.sqrt(0.125)) - 1e-12)
        self.assertLessEqual(twice_area.max()/2, self.options['max_area']*(1+1e-12))
        errors = triangle_errors(vertices, triangles, self.elevation, self.transform, self.srs,
                                 self.options['min_area'])
        self.assertGreater(len(errors), 0)
        maximum = errors.max()
        self.assertLessEqual(maximum, tolerance + 1e-7)
        return maximum

    # Route the shared fixture through a whole-domain solve or the selected seam merge, keeping
    # the terrain and refinement settings comparable.
    def mesh(self, options, bounds, on_disk):
        if bounds is None:
            return remesh_domain(self.domain, self.streams, options, self.work)
        tiles = [remesh_domain(tile_domain(self.domain, core, options['remesh_band_width']),
                               self.streams, options, self.work) for core in bounds]
        if on_disk:
            paths = [self.work/'left.npz', self.work/'right.npz']
            for path, (vertices, triangles) in zip(paths, tiles):
                np.savez(path, verts=vertices, tris=triangles)
            output = self.work/'merged.npz'
            merge_disk_cavity(paths[0], bounds[0], paths[1], bounds[1], self.domain,
                              self.streams, options, output)
            with np.load(output) as mesh:
                return mesh['verts'], mesh['tris']
        vertices, triangles, _ = merge_cavity(*tiles[0], bounds[0], *tiles[1], bounds[1],
                                             self.domain, self.streams, options, self.work)
        return vertices, triangles

    # Repeat the same quality checks from sub-metre to hundreds-of-metres tolerances. Check that
    # varying the tolerance actually changes the result.
    def run_sweep(self, bounds=None, on_disk=False):
        observed = {}
        for tolerance in TOLERANCES_M:
            with self.subTest(tolerance_m=tolerance):
                vertices, triangles = self.mesh(dict(self.options, max_tolerance=tolerance), bounds, on_disk)
                maximum = self.check_mesh(vertices, triangles, tolerance)
                observed[tolerance] = (len(triangles), maximum)
        if len(observed) == len(TOLERANCES_M):
            # Always over-refining also satisfies every upper bound. Verify
            # the fractional and coarse tolerances actually change the result.
            self.assertGreater(observed[0.75][1], 0.25)
            self.assertGreater(observed[500][1], 100)
            self.assertGreater(observed[0.25][0], observed[500][0])

    # Establish the single-domain result across tolerances without an artificial tile interface.
    def test_whole_domain(self):
        self.run_sweep()

    # Place the tile split across the domain and use the disk join so terrain quality is also
    # checked at the boundary of the retained interior.
    def test_vertical_disk_seam(self):
        x, y, xmax, ymax = self.domain.bounds
        self.run_sweep(([x, y, x+612, ymax], [x+612, y, xmax, ymax]), on_disk=True)

    # Rotate the artificial split and exercise the in-memory stitch so the sweep covers both
    # merge paths.
    def test_horizontal_in_memory_seam(self):
        x, y, xmax, ymax = self.domain.bounds
        self.run_sweep(([x, y, xmax, y+525], [x, y+525, xmax, ymax]))
