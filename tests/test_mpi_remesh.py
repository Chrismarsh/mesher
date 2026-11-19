import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import shutil
import subprocess
import sys

import numpy as np
from osgeo import gdal, osr
from shapely.geometry import box, Polygon, LineString, mapping, shape

from mesher.mesher_utls.mpi_remesh import (
    PatchError, merge_cavity, mesh_coverage, remesh_domain, validate_patch, _stitch,
    _mesher_arguments)
from mesher.mesh_topology import MeshTopologyError, audit_mesh

gdal.UseExceptions()


def grid(x0, x1, y0, y1, step=1):
    xs, ys = np.arange(x0, x1 + step / 2, step), np.arange(y0, y1 + step / 2, step)
    vertices = np.array([[x, y, 0] for y in ys for x in xs], dtype=float)
    tris = []
    for j in range(len(ys) - 1):
        for i in range(len(xs) - 1):
            a = j * len(xs) + i
            tris += [[a, a + 1, a + len(xs) + 1], [a, a + len(xs) + 1, a + len(xs)]]
    return vertices, np.array(tris)


def audit(vertices, triangles):
    from mesher.mesh_topology import build_neighbours
    return audit_mesh({'mesh': {'vertex': vertices.tolist(), 'elem': triangles.tolist(),
        'neigh': build_neighbours(triangles, len(vertices), vertices).tolist(),
        'nvertex': len(vertices), 'nelem': len(triangles)}})[0]


class CavityGeometryTests(unittest.TestCase):
    def test_directed_coverage_preserves_hole_and_disconnected_component(self):
        v, t = grid(0, 5, 0, 5)
        centers = v[t, :2].mean(axis=1)
        hole = ((centers > 1) & (centers < 3)).all(axis=1)
        t = t[~hole]
        w, u = grid(7, 8, 0, 1)
        coverage = mesh_coverage(np.vstack((v, w)), np.vstack((t, u + len(v))))
        expected = box(0, 0, 5, 5).difference(box(1, 1, 3, 3)).union(box(7, 0, 8, 1))
        self.assertEqual(coverage.symmetric_difference(expected).area, 0)

    def test_rejects_overlap_gap_and_degenerate_triangles(self):
        v, t = grid(0, 2, 0, 2)
        with self.assertRaisesRegex(PatchError, 'cover exactly'):
            validate_patch(v, t[:-1], box(0, 0, 2, 2))
        w = np.array([[.1, .1, 0], [.4, .1, 0], [.1, .4, 0]])
        with self.assertRaisesRegex(PatchError, 'overlap'):
            validate_patch(np.vstack((v, w)), np.vstack((t, [len(v), len(v)+1, len(v)+2])), box(0, 0, 2, 2))
        with self.assertRaisesRegex(MeshTopologyError, 'degenerate'):
            validate_patch(np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]), np.array([[0, 1, 2]]), box(0, 0, 2, 2))

    def test_split_retained_interface_is_rejected_even_with_exact_coverage(self):
        v = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
        t = np.array([[0, 1, 2], [0, 2, 3]])
        p = np.array([[1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [1, .5, 0]])
        q = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4]])
        validate_patch(p, q, box(1, 0, 2, 1))
        with self.assertRaisesRegex(PatchError, 'retained interface'):
            _stitch(v, t, p, q, box(0, 0, 2, 1))

    def test_unrepresented_physical_constraint_is_rejected(self):
        v, t = grid(0, 2, 0, 2)
        with self.assertRaisesRegex(PatchError, 'physical constraint'):
            validate_patch(v, t, box(0, 0, 2, 2), [LineString([(0, .3), (2, .3)])])

    def test_mixed_raster_weights_follow_cpp_reader_order(self):
        args = _mesher_arguments({'mesher_path': 'meshercxx', 'max_area': 8, 'min_area': .25,
            'errormetric': 'rmse', 'dem_path': 'dem.tif', 'max_tolerance': 1,
            'use_weights': True, 'topo_weight': .2, 'weight_threshold': .6,
            'parameter_files': {
                'landcover': {'filename': ['landcover.tif'], 'method': 'mode', 'tolerance': .8, 'weight': .5},
                'slope': {'filename': ['slope.tif'], 'method': 'mean', 'tolerance': 2, 'weight': .3}}})
        def values(flag): return [args[i+1] for i, arg in enumerate(args) if arg == flag]
        self.assertEqual(values('--raster'), ['dem.tif', 'slope.tif'])
        self.assertEqual(values('--category-raster'), ['landcover.tif'])
        self.assertEqual(values('--weight'), ['0.2', '0.3', '0.5'])

    def test_failed_merge_does_not_publish_candidate(self):
        from mesher.mesher_utls.MPI_merge_tiles import merge_tiles
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            inputs = []
            for name, bounds in (('a', [0, 0, 1, 1]), ('b', [1, 0, 2, 1])):
                np.savez(work / (name + '.npz'), verts=np.zeros((0, 3)), tris=np.zeros((0, 3), dtype=int))
                (work / (name + '.json')).write_text(json.dumps({'core_bbox': bounds}))
                inputs.append({'npz': str(work / (name + '.npz')), 'meta': str(work / (name + '.json'))})
            (work / 'joined.npz').write_bytes(b'previous mesh')
            (work / 'joined.json').write_bytes(b'previous metadata')
            with patch('mesher.mesher_utls.mpi_remesh.load_domain', return_value=box(0, 0, 2, 1)), \
                 patch('mesher.mesher_utls.mpi_remesh.merge_cavity', side_effect=PatchError('bad candidate')):
                with self.assertRaisesRegex(PatchError, 'bad candidate'):
                    merge_tiles({'tile_a': inputs[0], 'tile_b': inputs[1], 'merge_method': 'remesh',
                                 'outer_polygon_shp': 'unused', 'constraints': [], 'out_prefix': str(work / 'joined')})
            self.assertEqual((work / 'joined.npz').read_bytes(), b'previous mesh')
            self.assertEqual((work / 'joined.json').read_bytes(), b'previous metadata')


@unittest.skipUnless(os.environ.get('MESHER_TEST_EXE'), 'Set MESHER_TEST_EXE to the rebuilt meshercxx')
class CgalCavityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.dem = self.work / 'dem.tif'
        ds = gdal.GetDriverByName('GTiff').Create(str(self.dem), 160, 120, 1, gdal.GDT_Float64)
        ds.SetGeoTransform((-20, .5, 0, 40, 0, -.5))
        srs = osr.SpatialReference(); srs.ImportFromEPSG(32608)
        ds.SetProjection(srs.ExportToWkt()); ds.GetRasterBand(1).Fill(10); ds = None
        self.options = {'mesher_path': os.environ['MESHER_TEST_EXE'], 'max_area': 8,
            'min_area': .25, 'errormetric': 'rmse', 'max_tolerance': 1,
            'dem_path': str(self.dem), 'remesh_band_width': 2, 'remesh_max_attempts': 6,
            'remesh_max_vertices': 10000}

    def mesh(self, domain, constraints=()):
        return remesh_domain(domain, constraints, self.options, self.work)

    def test_backend_holes_components_and_crossing_constraints(self):
        domain = box(0, 0, 20, 20).difference(box(7, 7, 13, 13)).union(box(25, 0, 30, 5))
        lines = [LineString([(0, 5), (20, 5)]), LineString([(5, 0), (5, 20)])]
        v, t = self.mesh(domain, lines)
        validate_patch(v, t, domain, lines)
        self.assertTrue(audit(v, t)['valid'])
        self.assertEqual(audit(v, t)['components'], 2)

    def test_closed_physical_constraint_does_not_create_a_hole(self):
        domain = box(0, 0, 20, 20)
        lines = [LineString(box(5, 5, 15, 15).exterior.coords)]
        v, t = self.mesh(domain, lines)
        validate_patch(v, t, domain, lines)
        self.assertTrue(audit(v, t)['valid'])
        self.assertEqual(audit(v, t)['components'], 1)

    def test_rounded_exterior_centroid_cannot_start_runaway_refinement(self):
        # Reduced from the production boundary. A nearly collinear exterior
        # CDT face has a rounded centroid which OGR classifies as inside.
        # Geometry alone reproduced the runaway; no raster is required.
        fixture = Path(__file__).parent / 'fixtures' / 'raster_stair_boundary.geojson'
        domain = shape(json.loads(fixture.read_text()))
        constraints = self.work / 'constraints.json'
        constraints.write_text(json.dumps({'type': 'FeatureCollection', 'features': [
            {'type': 'Feature', 'properties': {}, 'geometry': mapping(LineString(domain.exterior.coords))}]}))
        output = self.work / 'geometry.json'
        result = subprocess.run([self.options['mesher_path'], '--area', '6250000', '--min-area', '2500',
            '--domain-file', str(fixture.resolve()), '--interior-plgs-file', str(constraints),
            '--patch-output', str(output), '--patch-max-vertices', '5000'],
            text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        mesh = json.loads(output.read_text())
        v, t = np.array(mesh['vertices']), np.array(mesh['triangles'])
        self.assertLess(len(v), 1000)
        validate_patch(v, t, domain)
        self.assertTrue(audit(v, t)['valid'])

    def test_overlapping_tiles_merge_with_no_seam_line_or_density_band(self):
        bounds_a, bounds_b = [0, 0, 20, 20], [20, 0, 40, 20]
        domain = box(0, 0, 40, 20)
        a, ta = self.mesh(box(0, 0, 23, 20))
        b, tb = self.mesh(box(17, 0, 40, 20))
        before_a, before_ta = a.copy(), ta.copy()
        v, t, stats = merge_cavity(a, ta, bounds_a, b, tb, bounds_b, domain, [], self.options, self.work)
        validate_patch(v, t, domain)
        self.assertTrue(audit(v, t)['valid'])
        self.assertEqual(audit(v, t)['hanging_nodes'], 0)
        np.testing.assert_array_equal(a, before_a); np.testing.assert_array_equal(ta, before_ta)
        self.assertGreater(stats['patch_triangles'], 0)
        # At most a few incidental points can lie on the now unconstrained seam.
        self.assertLess(np.count_nonzero((abs(v[:, 0] - 20) < 1e-9) & (v[:, 1] > 0) & (v[:, 1] < 20)), 8)
        serial_v, serial_t = self.mesh(domain)
        self.assertLess(len(t), 1.5 * len(serial_t))

    def test_retry_grows_cavity_and_exhaustion_does_not_modify_input(self):
        a, ta = grid(0, 20, 0, 20, 2)
        b, tb = grid(20, 40, 0, 20, 2)
        calls = []
        def retry_once(domain, constraints, options, directory, fixed):
            calls.append(domain.area)
            if len(calls) == 1:
                raise PatchError('Refinement changed a retained interface')
            return remesh_domain(domain, constraints, options, directory, fixed)
        v, t, stats = merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20],
            box(0, 0, 40, 20), [], self.options, self.work, retry_once)
        self.assertGreater(calls[1], calls[0]); self.assertGreater(stats['attempts'], 1)
        validate_patch(v, t, box(0, 0, 40, 20))
        original = ta.copy()
        def always_bad(*args): raise PatchError('invalid candidate')
        with self.assertRaisesRegex(PatchError, 'No valid cavity'):
            merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20],
                box(0, 0, 40, 20), [], dict(self.options, remesh_max_attempts=2), self.work, always_bad)
        np.testing.assert_array_equal(ta, original)

    def test_backend_budget_failure_is_explicit(self):
        with self.assertRaisesRegex(PatchError, 'budget.*vertices=4.*limit=4') as failure:
            remesh_domain(box(0, 0, 20, 20), [], dict(self.options, remesh_max_vertices=4), self.work)
        folders = list(self.work.glob('cavity-*'))
        self.assertEqual(len(folders), 1)
        self.assertIn(str(folders[0]), str(failure.exception))
        invocation = json.loads((folders[0] / 'invocation.json').read_text())
        self.assertEqual(invocation['domain_area'], 400)
        self.assertIn('refinement_steps=', (folders[0] / 'backend.log').read_text())
        self.assertTrue((folders[0] / 'domain.json').is_file())
        self.assertTrue((folders[0] / 'constraints.json').is_file())

    def test_backend_partial_and_disjoint_raster_coverage_is_bounded(self):
        # Source footprint is [-20, -20, 60, 40]. These include negative read
        # offsets and wholly disjoint windows which formerly underflowed in MEM.
        # This exercises sampling robustness, not valid DEM coverage for export.
        for bounds in [(-30, 0, -16, 20), (56, 0, 70, 20), (0, 36, 20, 50),
                       (0, -30, 20, -16), (-32, 0, -25, 20), (65, 0, 72, 20)]:
            with self.subTest(bounds=bounds):
                domain = box(*bounds)
                v, t = self.mesh(domain)
                validate_patch(v, t, domain)
                self.assertLess(len(v), 300)
        self.assertEqual(list(self.work.glob('cavity-*')), [])

    def test_four_tile_junction_with_crossing_stream_and_hole(self):
        domain = box(0, 0, 40, 20).difference(box(2, 2, 5, 5))
        constraints = [LineString([(0, 10), (40, 10)]), LineString([(20, 0), (20, 20)])]
        bounds = [[0, 0, 20, 10], [20, 0, 40, 10], [0, 10, 20, 20], [20, 10, 40, 20]]
        tiles = []
        for x0, y0, x1, y1 in bounds:
            tiles.append(self.mesh(domain.intersection(box(x0-3, y0-3, x1+3, y1+3)), constraints))
        rows = []
        for i in (0, 2):
            rows.append(merge_cavity(*tiles[i], bounds[i], *tiles[i+1], bounds[i+1],
                                    domain, constraints, self.options, self.work)[:2])
        v, t, _ = merge_cavity(*rows[0], [0, 0, 40, 10], *rows[1], [0, 10, 40, 20],
                               domain, constraints, self.options, self.work)
        validate_patch(v, t, domain, constraints)
        self.assertTrue(audit(v, t)['valid'])

    def test_empty_tile_keeps_domain_and_ownership_boxes(self):
        a, ta = self.mesh(box(1, 1, 10, 19))
        v, t, _ = merge_cavity(a, ta, [0, 0, 20, 20], np.empty((0, 3)), np.empty((0, 3), dtype=int),
            [20, 0, 40, 20], box(1, 1, 10, 19), [], self.options, self.work)
        validate_patch(v, t, box(1, 1, 10, 19))

    def test_rugged_dem_is_refined_by_original_error_criteria(self):
        ds = gdal.Open(str(self.dem), gdal.GA_Update)
        x = -20 + (np.arange(160) + .5) * .5
        y = 40 - (np.arange(120) + .5) * .5
        xx, yy = np.meshgrid(x, y)
        ds.GetRasterBand(1).WriteArray(10 + 30 * np.exp(-((xx-20)**2 + (yy-10)**2) / 20))
        ds = None
        domain = box(0, 0, 40, 20)
        a, ta = self.mesh(box(0, 0, 23, 20)); b, tb = self.mesh(box(17, 0, 40, 20))
        v, t, _ = merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20],
                               domain, [], self.options, self.work)
        validate_patch(v, t, domain)
        self.assertTrue(audit(v, t)['valid'])

    def test_shifted_seam_sloping_dem_and_diagonal_constraint(self):
        ds = gdal.Open(str(self.dem), gdal.GA_Update)
        xx, yy = np.meshgrid(-20 + (np.arange(160) + .5) * .5, 40 - (np.arange(120) + .5) * .5)
        ds.GetRasterBand(1).WriteArray(10 + .5 * xx + .3 * yy); ds = None
        domain = box(0, 0, 40, 20)
        constraints = [LineString([(5, 0), (35, 20)])]
        a, ta = self.mesh(box(0, 0, 20, 20), constraints)
        b, tb = self.mesh(box(14, 0, 40, 20), constraints)
        v, t, _ = merge_cavity(a, ta, [0, 0, 17, 20], b, tb, [17, 0, 40, 20],
                               domain, constraints, self.options, self.work)
        validate_patch(v, t, domain, constraints)
        self.assertTrue(audit(v, t)['valid'])

    def test_projected_coordinate_magnitudes_and_diagonal_constraint(self):
        from shapely.affinity import translate
        x, y = 500000, 7600000
        ds = gdal.Open(str(self.dem), gdal.GA_Update)
        ds.SetGeoTransform((x-20, .5, 0, y+40, 0, -.5)); ds = None
        domain = translate(box(0, 0, 40, 20), x, y)
        constraints = [translate(LineString([(5, 0), (35, 20)]), x, y)]
        a, ta = self.mesh(translate(box(0, 0, 20, 20), x, y), constraints)
        b, tb = self.mesh(translate(box(14, 0, 40, 20), x, y), constraints)
        v, t, _ = merge_cavity(a, ta, [x, y, x+17, y+20], b, tb, [x+17, y, x+40, y+20],
                               domain, constraints, self.options, self.work)
        validate_patch(v, t, domain, constraints)
        self.assertTrue(audit(v, t)['valid'])

    @unittest.skipUnless(shutil.which('mpirun'), 'OpenMPI is required')
    def test_real_four_rank_pipeline_and_export(self):
        from mesher.mesher_cli import run_mpi_meshing, build_mesh_from_arrays
        domain_file = self.work / 'domain.geojson'
        domain_file.write_text(json.dumps({'type': 'Polygon', 'coordinates':
            [[[0, 0], [40, 0], [40, 20], [0, 20], [0, 0]]]}))
        shutil.copyfile(self.dem, self.work / 'case_projected.tif')
        cwd = os.getcwd()
        try:
            os.chdir(self.work)
            final = run_mpi_meshing(str(self.work) + '/', 'case', 0, 0, 40, 20,
                '/usr/bin', self.options['mesher_path'], str(domain_file), {}, {}, {},
                8, .25, 1, 'rmse', 0, False, 1, 0, False,
                f'mpirun --oversubscribe -n 4 {sys.executable}', 4,
                remesh_options={'merge_method': 'remesh', 'remesh_band_width': 2,
                                'remesh_max_attempts': 6, 'remesh_max_vertices': 10000})
        finally:
            os.chdir(cwd)
        with np.load(final) as result:
            v, t = result['verts'], result['tris']
        validate_patch(v, t, box(0, 0, 40, 20))
        ds = gdal.Open(str(self.dem)); srs = osr.SpatialReference(); srs.ImportFromEPSG(32608)
        exported = build_mesh_from_arrays(v, t, ds, ds.GetRasterBand(1), srs, False, False)
        self.assertTrue(audit_mesh(exported)[0]['valid'])
        self.assertEqual(exported['mesh']['nelem'], len(t))
