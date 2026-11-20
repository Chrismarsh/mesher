import json
import io
from contextlib import redirect_stdout
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
from shapely import covers, polygons
from shapely.geometry import box, Polygon, LineString, mapping, shape
from shapely.affinity import rotate, translate, scale

from mesher.mesher_utls.mpi_remesh import (
    PatchError, PatchQualityError, merge_cavity, mesh_coverage, remesh_domain,
    validate_patch, _stitch, _mesher_arguments, _quality_interface_faces,
    InterfaceChanged, _grow_cavity, tile_domain)
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


def buffered_slanted_domain():
    # The physical rectangle has right angles; rectangular clipping creates
    # the acute 20.1093-degree corner reported for the five-vertex tile input.
    rectangle = translate(rotate(box(-50, -2.5, 50, 2.5), 20.1093, origin=(0, 0)), -22, -17)
    rectangle = scale(rectangle, xfact=1000, yfact=1000, origin=(0, 0))
    return rectangle.buffer(200, join_style=2).simplify(10, preserve_topology=True).buffer(-200)


class CavityGeometryTests(unittest.TestCase):
    def test_local_cut_extensions_keep_owned_area_holes_and_physical_boundary(self):
        for x, y in ((0, 0), (-1934381.4845238777, 3429022.4631309444)):
            with self.subTest(origin=(x, y)):
                domain = translate(buffered_slanted_domain(), x, y)
                core = box(x, y, x+20000, y+20000)
                centre = domain.intersection(core).buffer(-50).representative_point()
                hole = box(centre.x-25, centre.y-25, centre.x+25, centre.y+25)
                self.assertTrue(domain.intersection(core).covers(hole))
                domain = domain.difference(hole)
                raw = domain.intersection(box(x-5000, y-5000, x+25000, y+25000))
                result = tile_domain(domain, list(core.bounds), 5000)
                self.assertGreater(result.area, raw.area)
                self.assertLess(result.area-raw.area, 5000**2)
                self.assertEqual(domain.intersection(core).difference(result).area, 0)
                self.assertAlmostEqual(result.difference(domain).area, 0, places=5)
                self.assertEqual(result.intersection(hole).area, 0)
                self.assertEqual(result.difference(raw).intersection(core).area, 0)

    def test_halo_point_or_line_contact_does_not_own_core_area(self):
        for domain in (box(-2, 2, -1, 3), box(-2, 2, 0, 3), box(-2, -2, 0, 0)):
            with self.subTest(domain=domain.wkt):
                self.assertTrue(tile_domain(domain, [0, 0, 10, 10], 3).is_empty)

    def test_tile_keeps_owned_halo_and_holes_but_omits_unowned_components(self):
        owned = box(-1, -1, 4, 4).difference(box(1, 1, 2, 2))
        outside = box(11, 1, 13, 2)
        touching = box(-1, 6, 0, 7)
        domain = owned.union(outside).union(touching)
        result = tile_domain(domain, [0, 0, 10, 10], 2)
        self.assertEqual(result.symmetric_difference(owned).area, 0)
        self.assertEqual(len(result.interiors), 1)
        # Dropped components remain in the unchanged global domain.
        self.assertTrue(domain.covers(outside.union(touching)))

    def test_arbitrarily_thin_owned_component_keeps_its_full_halo(self):
        main = box(1, 1, 3, 3)
        thin = box(-1, 6, 1e-9, 7)
        domain = main.union(thin)
        result = tile_domain(domain, [0, 0, 10, 10], 2)
        self.assertGreater(thin.intersection(box(0, 0, 10, 10)).area, 0)
        self.assertEqual(result.symmetric_difference(domain).area, 0)

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
        for x, y in ((0, 0), (-1934381.4845238777, 3429022.4631309444)):
            with self.subTest(origin=(x, y)):
                v = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]) + [x, y, 0]
                t = np.array([[0, 1, 2], [0, 2, 3]])
                p = np.array([[1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [1, .5, 0]]) + [x, y, 0]
                q = np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4]])
                validate_patch(p, q, box(x+1, y, x+2, y+1))
                with self.assertRaisesRegex(PatchError, 'retained interface'):
                    _stitch(v, t, p, q, box(x, y, x+2, y+1))

    def test_boundary_roundoff_allowance_does_not_accept_a_real_gap(self):
        x, y = -1934381.4845238777, 3429022.4631309444
        v = np.array([[x, y, 0], [x+1, y, 0], [x+1, y+1-1e-5, 0], [x, y+1-1e-5, 0]])
        t = np.array([[0, 1, 2], [0, 2, 3]])
        with self.assertRaises(InterfaceChanged):
            _stitch(v, t, np.empty((0, 3)), np.empty((0, 3), dtype=int), box(x, y, x+1, y+1))

    def test_local_growth_does_not_cross_existing_cavities(self):
        # A face chain separated by an already removed face. Even a deep local
        # expansion must stop at that gap and preserve the other component.
        neighbours = np.full((20, 3), -1, dtype=int)
        neighbours[:-1, 0] = np.arange(1, 20)
        neighbours[1:, 1] = np.arange(19)
        keep = np.ones(20, dtype=bool)
        keep[10] = False
        self.assertEqual(_grow_cavity(keep, neighbours, [0], 100), 10)
        np.testing.assert_array_equal(np.flatnonzero(keep), np.arange(11, 20))

    def test_repeated_interface_splits_expand_faster_than_one_ring(self):
        # Each valid finer patch bisects the retained grid's edges. A one-ring
        # retry keeps moving the mismatch and exhausts eight attempts. Growing
        # local adjacency depth resolves it while retaining a separate island.
        a, ta = grid(0, 50, 0, 4)
        island_v, island_t = grid(0, 4, 10, 14)
        ta = np.vstack((ta, island_t + len(a)))
        a = np.vstack((a, island_v))
        b, tb = grid(50, 100, 0, 4)
        original_a, original_ta = a.copy(), ta.copy()
        domain = box(0, 0, 100, 4).union(box(0, 10, 4, 14))
        fine_v, fine_t = grid(0, 100, 0, 4, .5)
        fine_polygons = polygons(fine_v[fine_t, :2])
        for limit in (4, 10):
            with self.subTest(attempt_limit=limit):
                calls, interfaces = [], []
                def finer_patch(cavity, constraints, options, directory, fixed):
                    calls.append(cavity.area)
                    interfaces.append(len(fixed))
                    selected = fine_t[covers(cavity, fine_polygons)]
                    used, inverse = np.unique(selected, return_inverse=True)
                    return fine_v[used], inverse.reshape((-1, 3))
                v, t, stats = merge_cavity(a, ta, [0, 0, 50, 14], b, tb, [50, 0, 100, 14], domain, [],
                    {'max_area': .5, 'remesh_band_width': 1, 'remesh_max_attempts': limit}, '.', finer_patch)
                self.assertGreater(stats['attempts'], 2)
                self.assertLessEqual(stats['attempts'], limit)
                self.assertEqual(stats['band_width'], 1)
                self.assertEqual(stats['retained_triangles'], len(island_t))
                self.assertEqual(stats['interface_recovery'], limit == 4)
                self.assertEqual(interfaces[-1], 0)
                self.assertTrue(all(second > first for first, second in zip(calls, calls[1:])))
                # Preserve every original triangle of the complete, separate
                # island even when the last attempt removes other interfaces.
                np.testing.assert_array_equal(v[t[:len(island_t)]], island_v[island_t])
                np.testing.assert_array_equal(a, original_a)
                np.testing.assert_array_equal(ta, original_ta)
                validate_patch(v, t, domain)
                report = audit(v, t)
                self.assertTrue(report['valid'], report)
                self.assertEqual(report['components'], 2)

    def test_interface_recovery_cannot_skip_validation_or_vertex_budget(self):
        a, ta = grid(0, 20, 0, 20, 2)
        b, tb = grid(20, 40, 0, 20, 2)
        before_a, before_b = ta.copy(), tb.copy()
        for failure in ('budget', 'gap', 'generic_gap'):
            with self.subTest(failure=failure):
                calls = []
                def candidate(cavity, constraints, options, directory, fixed):
                    calls.append(len(fixed))
                    if fixed:
                        if failure == 'generic_gap':
                            raise PatchError('Invalid candidate coverage')
                        raise InterfaceChanged([0])
                    if failure == 'budget':
                        raise PatchError('Patch vertex budget exceeded: vertices=100, limit=100')
                    return grid(0, 1, 0, 1)  # Invalid coverage must still fail.
                with self.assertRaises(PatchError) as error:
                    merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20], box(0, 0, 40, 20), [],
                        {'max_area': 8, 'remesh_band_width': 2, 'remesh_max_attempts': 2}, '.', candidate)
                self.assertEqual(len(calls), 2)
                self.assertEqual(calls[-1], 0)
                self.assertIn('after 2 attempts', str(error.exception))
                self.assertIn('vertex budget exceeded' if failure == 'budget' else 'cover exactly', str(error.exception))
                np.testing.assert_array_equal(ta, before_a)
                np.testing.assert_array_equal(tb, before_b)

    def test_last_attempt_does_not_grow_an_untested_cavity(self):
        a, ta = grid(0, 20, 0, 20, 2)
        b, tb = grid(20, 40, 0, 20, 2)
        def changed(*args): raise InterfaceChanged([0])
        log = io.StringIO()
        with redirect_stdout(log), patch('mesher.mesher_utls.mpi_remesh._grow_cavity') as grow:
            with self.assertRaisesRegex(PatchError, 'after 1 attempts.*stage0_r4_c0'):
                merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20], box(0, 0, 40, 20), [],
                    {'max_area': 8, 'remesh_band_width': 2, 'remesh_max_attempts': 1,
                     'out_prefix': '/example/stage0_r4_c0'}, '.', changed)
        grow.assert_not_called()
        self.assertNotIn('Growing cavity', log.getvalue())
        self.assertIn('[stage0_r4_c0]', log.getvalue())

    def test_unrepresented_physical_constraint_is_rejected(self):
        v, t = grid(0, 2, 0, 2)
        with self.assertRaisesRegex(PatchError, 'physical constraint'):
            validate_patch(v, t, box(0, 0, 2, 2), [LineString([(0, .3), (2, .3)])])

    def test_quality_failure_selects_only_touching_interface_owners_at_large_coordinates(self):
        # From the production cavity: a refined point is about 2e-10 m off
        # the sloping interface after rounding. Selection may tolerate that;
        # stitching must still require exact, unsplit edges.
        tip = (-1983662.301079662, 3508355.7468695333)
        edges = [LineString([(-1982502.0452924795, 3507007.5602487084), tip]),
                 LineString([tip, (-1984237.619732993, 3506687.2707918254)]),
                 LineString([(-1990000, 3500000), (-1991000, 3501000)])]
        failures = [{'vertices': [(-1984335.8080216227, 3508466.0784428404),
                                 (-1983949.9604063276, 3507521.5088306796),
                                 (-1984004.6785210543, 3508466.0784428404)]}]
        domain = box(-2010827, 3434022, -1841287, 3508467)
        np.testing.assert_array_equal(_quality_interface_faces(edges, [3, 7, 9], failures, domain), [7])
        failures.append({'vertices': [(-1984004.6785210543, 3508466.0784428404), tip,
                                     (-1983839.1137707701, 3508466.0784428404)]})
        np.testing.assert_array_equal(_quality_interface_faces(edges, [3, 7, 9], failures, domain), [3, 7])
        self.assertEqual(len(_quality_interface_faces([], [], failures, domain)), 0)

    def test_vertex_budget_failure_does_not_retry_a_larger_cavity(self):
        a, ta = grid(0, 20, 0, 20, 2)
        b, tb = grid(20, 40, 0, 20, 2)
        with patch('mesher.mesher_utls.mpi_remesh.remesh_domain',
                   side_effect=PatchError('Patch vertex budget exceeded: vertices=500000, limit=500000')) as backend:
            with self.assertRaisesRegex(PatchError, 'after 1 attempts.*vertex budget exceeded'):
                merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20], box(0, 0, 40, 20), [],
                             {'max_area': 8, 'remesh_band_width': 2}, '.', backend)
        self.assertEqual(backend.call_count, 1)

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

    def test_buffered_halo_only_triangle_does_not_abort_empty_core_tile(self):
        from mesher.mesher_utls.MPI_mesh_tiles import mesh_tile
        angle = np.deg2rad(19.2648)
        height = np.sqrt(2*3338810*np.sin(angle)*np.cos(angle))
        p = np.array([-5000+height, 10000])
        u = 6000*np.array([-np.cos(angle), np.sin(angle)])
        w = 10000*np.array([-np.sin(angle), -np.cos(angle)])
        raw = Polygon([p, p+u, p+u+w, p+w]).buffer(200, join_style=2)
        domain = raw.simplify(10, preserve_topology=True).buffer(-200)
        halo_fragment = domain.intersection(box(-5000, -5000, 25000, 25000))
        self.assertEqual(len(halo_fragment.exterior.coords)-1, 3)
        self.assertAlmostEqual(halo_fragment.area, 3338810)
        self.assertEqual(domain.intersection(box(0, 0, 20000, 20000)).area, 0)
        options = dict(self.options, max_area=6250000, min_area=2500, max_tolerance=-1,
                       remesh_band_width=5000, remesh_max_vertices=500000)
        with self.assertRaises(PatchQualityError) as failure:
            remesh_domain(halo_fragment, [], options, self.work)
        self.assertIn('vertices=3', str(failure.exception))
        self.assertTrue(any(f['reason'] == 'angle' and abs(f['min_angle_degrees']-19.2648) < 1e-4
                            for f in failure.exception.failures))
        domain_file = self.work / 'buffered.geojson'
        domain_file.write_text(json.dumps(mapping(domain)))
        args = dict(options, merge_method='remesh', tile_bbox=[0, 0, 20000, 20000],
                    band_width=5000, outer_polygon_shp=str(domain_file), constraints=[],
                    base_dir=str(self.work)+'/', tile_prefix='unowned')
        mesh_tile(args)
        metadata = json.loads((self.work / 'unowned.json').read_text())
        self.assertTrue(metadata['empty_tile'])
        self.assertEqual(metadata['core_bbox'], args['tile_bbox'])
        with np.load(self.work / 'unowned.npz') as empty:
            self.assertEqual(empty['verts'].shape, (0, 3))
            self.assertEqual(empty['tris'].shape, (0, 3))
            # The neighbour owns the real domain. Omitting redundant halo work
            # must not omit any of that domain from the eventual merged mesh.
            owner_args = dict(args, tile_bbox=[-20000, 0, 0, 20000], tile_prefix='owner')
            mesh_tile(owner_args)
            with np.load(self.work / 'owner.npz') as owner:
                v, t, _ = merge_cavity(owner['verts'], owner['tris'], owner_args['tile_bbox'],
                    empty['verts'], empty['tris'], args['tile_bbox'], domain, [], options, self.work)
        validate_patch(v, t, domain)
        self.assertTrue(audit(v, t)['valid'])

    def test_owned_acute_physical_corner_still_fails_tile_quality(self):
        from mesher.mesher_utls.MPI_mesh_tiles import mesh_tile
        domain_file = self.work / 'acute.geojson'
        domain_file.write_text(json.dumps(mapping(Polygon([(1, 1), (15, 1), (1, 3)]))))
        args = dict(self.options, merge_method='remesh', tile_bbox=[0, 0, 20, 20], band_width=2,
                    outer_polygon_shp=str(domain_file), constraints=[], base_dir=str(self.work)+'/',
                    tile_prefix='owned', max_area=100, min_area=100, max_tolerance=-1)
        with self.assertRaises(PatchQualityError):
            mesh_tile(args)
        self.assertFalse((self.work / 'owned.npz').exists())
        self.assertFalse((self.work / 'owned.json').exists())

    def test_five_vertex_owned_tile_extends_only_its_artificial_corner(self):
        from mesher.mesher_utls.MPI_mesh_tiles import mesh_tile
        options = dict(self.options, max_area=6250000, min_area=2500, max_tolerance=-1,
                       remesh_band_width=5000, remesh_max_vertices=500000)
        for x, y in ((0, 0), (-1934381.4845238777, 3429022.4631309444)):
            with self.subTest(origin=(x, y)):
                domain = translate(buffered_slanted_domain(), x, y)
                bounds = [x, y, x+20000, y+20000]
                raw = domain.intersection(box(x-5000, y-5000, x+25000, y+25000))
                self.assertEqual(len(raw.exterior.coords)-1, 5)
                self.assertGreater(domain.intersection(box(*bounds)).area, 0)
                with self.assertRaises(PatchQualityError) as failure:
                    remesh_domain(raw, [], options, self.work)
                self.assertIn('vertices=5', str(failure.exception))
                self.assertTrue(any(f['reason'] == 'angle' and abs(f['min_angle_degrees']-20.1093) < 1e-4
                                    for f in failure.exception.failures))
                domain_file = self.work / 'buffered.geojson'
                domain_file.write_text(json.dumps(mapping(domain)))
                args = dict(options, merge_method='remesh', tile_bbox=bounds, band_width=5000,
                            outer_polygon_shp=str(domain_file), constraints=[], base_dir=str(self.work)+'/',
                            tile_prefix='owned_cut')
                mesh_tile(args)
                with np.load(self.work / 'owned_cut.npz') as mesh:
                    v, t = mesh['verts'], mesh['tris']
                    validate_patch(v, t, tile_domain(domain, bounds, 5000))
                    self.assertTrue(audit(v, t)['valid'])
                    covered_core = mesh_coverage(v, t).intersection(box(*bounds))
                    self.assertLess(covered_core.symmetric_difference(domain.intersection(box(*bounds))).area, 1e-5)
                self.assertFalse(json.loads((self.work / 'owned_cut.json').read_text())['empty_tile'])

    def test_intermediate_cut_halos_disappear_in_complete_merge(self):
        from mesher.mesher_utls.MPI_mesh_tiles import mesh_tile
        from mesher.mesher_utls.MPI_merge_tiles import merge_tiles
        domain = buffered_slanted_domain()
        domain_file = self.work / 'whole.geojson'
        domain_file.write_text(json.dumps(mapping(domain)))
        options = dict(self.options, max_area=6250000, min_area=2500, max_tolerance=-1,
                       remesh_band_width=5000, remesh_max_vertices=500000, merge_method='remesh',
                       outer_polygon_shp=str(domain_file), constraints=[])
        bounds = [[-80000, -40000, 0, 0], [0, -40000, 40000, 0],
                  [-80000, 0, 0, 20000], [0, 0, 40000, 20000]]
        def files(name):
            return {'npz': str(self.work / (name+'.npz')), 'meta': str(self.work / (name+'.json'))}
        for i, core in enumerate(bounds):
            mesh_tile(dict(options, tile_bbox=core, band_width=5000,
                           base_dir=str(self.work)+'/', tile_prefix=f'tile{i}'))
        group_halo_area = []
        for row in range(2):
            name = f'row{row}'
            merge_tiles(dict(options, tile_a=files(f'tile{2*row}'), tile_b=files(f'tile{2*row+1}'),
                             out_prefix=str(self.work / name)))
            metadata = json.loads((self.work / (name+'.json')).read_text())
            core = box(*metadata['core_bbox'])
            group_halo_area.append(metadata['remesh_report']['boundary_halo_area'])
            with np.load(self.work / (name+'.npz')) as mesh:
                covered = mesh_coverage(mesh['verts'], mesh['tris'])
                self.assertLess(covered.intersection(core).symmetric_difference(domain.intersection(core)).area, 1e-5)
                self.assertTrue(audit(mesh['verts'], mesh['tris'])['valid'])
        self.assertGreater(max(group_halo_area), 0)
        merge_tiles(dict(options, tile_a=files('row0'), tile_b=files('row1'), out_prefix=str(self.work / 'final')))
        with np.load(self.work / 'final.npz') as mesh:
            validate_patch(mesh['verts'], mesh['tris'], domain)
            report = audit(mesh['verts'], mesh['tris'])
            self.assertTrue(report['valid'], report)
            self.assertEqual(report['components'], 1)
        metadata = json.loads((self.work / 'final.json').read_text())['remesh_report']
        self.assertAlmostEqual(metadata['boundary_halo_area'], 0)
        self.assertGreater(metadata['original_triangles_retained'], 0)

    def test_batched_subdivision_multiple_sides_at_large_coordinates(self):
        from mesher.mesher_utls.mpi_remesh import subdivide_parents
        records = []
        for i, (x, y) in enumerate(((0, 0), (-1934381.4845238777, 3429022.4631309444))):
            parent = np.array([[0, 0], [6, 0], [3, 3*np.sqrt(3)]]) + [x, y]
            boundary = [((1-f)*a + f*b).tolist() for a, b in zip(parent, np.roll(parent, -1, axis=0))
                        for f in (0, 1/3, 2/3)]
            records.append({'id': i, 'parent': parent.tolist(), 'boundary': boundary})
        result = subdivide_parents(records, dict(self.options, max_area=100, min_area=100), self.work)
        self.assertTrue(all(r['accepted'] for r in result['repairs']), result)
        self.assertEqual(result['stats']['raster_evaluations'], 0)
        self.assertEqual(len(list(self.work.glob('repair-*'))), 0)

    def test_small_children_skip_rasters_but_bad_shape_is_rejected(self):
        from mesher.mesher_utls.mpi_remesh import subdivide_parents
        records = [
            {'id': 7, 'parent': [[0, 0], [2, 0], [0, 2]],
             'boundary': [[0, 0], [2, 0], [1, 1], [0, 2]]},
            {'id': 9, 'parent': [[0, 0], [2, 0], [0, 2]],
             'boundary': [[0, 0], [.01, 0], [2, 0], [0, 2]]}]
        # tol=0 rejects even a perfectly flat raster above min_area.
        result = subdivide_parents(records, dict(self.options, min_area=2, max_tolerance=0), self.work)
        self.assertEqual([r['accepted'] for r in result['repairs']], [True, False])
        self.assertEqual(result['stats']['raster_evaluations'], 0)
        above = subdivide_parents(records[:1], dict(self.options, min_area=.25, max_tolerance=0), self.work)
        self.assertFalse(above['repairs'][0]['accepted'])
        self.assertGreater(above['stats']['raster_evaluations'], 0)

    def test_patch_interiors_survive_repair_and_later_failure_is_transactional(self):
        ds = gdal.Open(str(self.dem), gdal.GA_Update)
        xx, yy = np.meshgrid(-20+(np.arange(160)+.5)*.5, 40-(np.arange(120)+.5)*.5)
        ds.GetRasterBand(1).WriteArray(10+30*np.exp(-((xx-20)**2+(yy-10)**2)/20)); ds = None
        domain = box(0, 0, 40, 20)
        a, ta = self.mesh(box(0, 0, 23, 20)); b, tb = self.mesh(box(17, 0, 40, 20))
        originals = [x.copy() for x in (a, ta, b, tb)]
        calls, results = [], []
        def capture(cavity, constraints, options, directory, fixed):
            calls.append(cavity)
            result = remesh_domain(cavity, constraints, options, directory, fixed)
            results.append(result)
            return result
        v, t, stats = merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20],
                                  domain, [], self.options, self.work, capture)
        self.assertGreater(stats['generated_triangles_reused'], 0)
        self.assertGreater(stats['local_repairs'], 0)
        self.assertFalse(stats['interface_recovery'])
        self.assertLess(calls[-1].area, calls[0].area)
        def keys(vertices, triangles):
            return {tuple(sorted(map(tuple, face))) for face in vertices[triangles]}
        final = keys(v, t)
        first_v, first_t = results[0]
        untouched = [tri for tri in first_t if all(Polygon(first_v[tri, :2]).disjoint(region)
                                                  for region in calls[1:])]
        self.assertGreater(len(untouched), 0)
        self.assertTrue(keys(first_v, np.asarray(untouched)).issubset(final))
        original_untouched = [tri for tri in ta if box(0, 0, 20, 20).covers(Polygon(a[tri, :2]))
                              and all(Polygon(a[tri, :2]).disjoint(region) for region in calls)]
        self.assertTrue(keys(a, np.asarray(original_untouched)).issubset(final))
        validate_patch(v, t, domain)
        self.assertTrue(audit(v, t)['valid'])
        # Fail after promotion of the first valid patch. No input array changes.
        invocations = []
        def fail_after_preservation(*args):
            invocations.append(1)
            if len(invocations) > 1:
                raise PatchError('Patch vertex budget exceeded in later repair')
            return remesh_domain(*args)
        with self.assertRaisesRegex(PatchError, 'after 2 attempts.*budget'):
            merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20], domain, [],
                         self.options, self.work, fail_after_preservation)
        for actual, expected in zip((a, ta, b, tb), originals):
            np.testing.assert_array_equal(actual, expected)

    def test_backend_signal_stops_retry_and_keeps_diagnostics(self):
        from mesher.mesher_utls.mpi_remesh import BackendTerminated
        executable = self.work / 'terminated-backend'
        executable.write_text('#!/bin/sh\nkill -TERM "$$"\n')
        executable.chmod(0o755)
        a, ta = grid(0, 20, 0, 20, 2); b, tb = grid(20, 40, 0, 20, 2)
        with self.assertRaisesRegex(PatchError, 'after 1 attempts.*signal 15'):
            merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20], box(0, 0, 40, 20), [],
                         dict(self.options, mesher_path=str(executable)), self.work)
        self.assertEqual(len(list(self.work.glob('cavity-*'))), 1)

    def test_bad_subdivision_output_is_rejected_before_promotion(self):
        from mesher.mesher_utls.mpi_remesh import subdivide_parents
        record = {'id': 0, 'parent': [[0, 0], [2, 0], [0, 2]],
                  'boundary': [[0, 0], [2, 0], [1, 1], [0, 2]]}
        def wrong_output(command, work, options):
            (work / 'mesh.json').write_text(json.dumps({'repairs': [
                {'id': 0, 'accepted': True, 'triangles': [[[0, 0], [2, 0], [0, 1]]]}]}))
            return 0, .01
        with patch('mesher.mesher_utls.mpi_remesh._run_backend', side_effect=wrong_output):
            with self.assertRaisesRegex(PatchError, 'cover exactly'):
                subdivide_parents([record], self.options, self.work)
        self.assertEqual(len(list(self.work.glob('repair-*'))), 1)

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

    def notch_options(self):
        ds = gdal.Open(str(self.dem), gdal.GA_Update)
        ds.SetGeoTransform((0, 50, 0, 6000, 0, -50)); ds = None
        return dict(self.options, max_area=6250000, min_area=2500, remesh_band_width=1000)

    def test_quality_report_preserves_locations_and_rejects_bad_mesh(self):
        options = self.notch_options()
        notch = Polygon([(1000, 0), (3000, 0), (2000, 3450)])
        domain = box(0, 0, 8000, 3500).difference(notch)
        with self.assertRaises(PatchQualityError) as failure:
            remesh_domain(domain, [], options, self.work)
        folders = list(self.work.glob('cavity-*'))
        self.assertEqual(len(folders), 1)
        self.assertIn(str(folders[0]), str(failure.exception))
        self.assertFalse((folders[0] / 'mesh.json').exists())
        report = json.loads((folders[0] / 'mesh.json.quality.json').read_text())
        self.assertEqual(failure.exception.failures, report['faces'])
        self.assertGreater(report['failed_faces'], 0)
        self.assertEqual(report['failed_faces'], report['angle_failures'])
        self.assertEqual(report['area_failures'] + report['raster_failures'], 0)
        for face in report['faces']:
            self.assertEqual(face['reason'], 'angle')
            self.assertLess(face['min_angle_degrees'], np.degrees(np.arcsin(np.sqrt(.125))))
            self.assertAlmostEqual(Polygon(face['vertices']).area, face['area'])
            self.assertEqual(len(face['constrained_edges']), 3)
        self.assertIn('reason=angle', (folders[0] / 'backend.log').read_text())

    def test_quality_retry_absorbs_protrusion_without_widening_the_entire_band(self):
        options = self.notch_options()
        a = np.array([[1000, 0, 0], [3000, 0, 0], [2000, 3450, 0]], dtype=float)
        b = np.array([[6000, 0, 0], [8000, 0, 0], [8000, 2000, 0]], dtype=float)
        tris = np.array([[0, 1, 2]])
        domain = box(0, 0, 8000, 3500)
        cavities = []
        def record_cavity(cavity, constraints, args, directory, fixed):
            cavities.append(cavity)
            return remesh_domain(cavity, constraints, args, directory, fixed)
        v, t, stats = merge_cavity(a, tris, [0, 0, 4000, 3500], b, tris, [4000, 0, 8000, 3500],
                                  domain, [], options, self.work, record_cavity)
        self.assertEqual(stats['attempts'], 2)
        self.assertEqual(stats['band_width'], 1000)
        self.assertEqual(stats['retained_triangles'], 1)
        self.assertEqual(cavities[1].difference(cavities[0]).symmetric_difference(Polygon(a[:, :2])).area, 0)
        self.assertEqual(set(map(tuple, v[t[0]])), set(map(tuple, b)))
        validate_patch(v, t, domain)
        report = audit(v, t)
        self.assertTrue(report['valid'], report)
        self.assertEqual(report['hanging_nodes'], 0)

    def test_quality_retry_cannot_remove_a_physical_notch(self):
        options = self.notch_options()
        domain = box(0, 0, 8000, 3500).difference(Polygon([(1000, 0), (3000, 0), (2000, 3450)]))
        empty_v, empty_t = np.empty((0, 3)), np.empty((0, 3), dtype=int)
        with self.assertRaisesRegex(PatchError, 'No valid cavity replacement after 1 attempts'):
            merge_cavity(empty_v, empty_t, [0, 0, 4000, 3500], empty_v, empty_t, [4000, 0, 8000, 3500],
                         domain, [], options, self.work)
        folder, = self.work.glob('cavity-*')
        self.assertFalse((folder / 'mesh.json').exists())
        supplied = shape(json.loads((folder / 'domain.json').read_text()))
        self.assertEqual(supplied.symmetric_difference(domain).area, 0)

    def test_final_attempt_removes_interfaces_before_real_cgal_refinement(self):
        a, ta = grid(0, 20, 0, 20, 2)
        b, tb = grid(20, 40, 0, 20, 2)
        domain = box(0, 0, 40, 20).difference(box(26, 4, 30, 8))
        tb = tb[covers(domain, polygons(b[tb, :2]))]
        line = LineString([(0, 0), (40, 20)])
        calls = []
        def split_interface(cavity, constraints, options, directory, fixed):
            calls.append(len(fixed))
            forced = list(fixed)
            if forced:
                edge = forced[0]
                # A valid refined patch which always subdivides a retained
                # edge reproduces the endless moving-interface failure.
                forced[0] = LineString([edge.coords[0], edge.interpolate(.5, normalized=True).coords[0], edge.coords[-1]])
            return remesh_domain(cavity, constraints, options, directory, forced)
        v, t, stats = merge_cavity(a, ta, [0, 0, 20, 20], b, tb, [20, 0, 40, 20], domain, [line],
            dict(self.options, remesh_max_attempts=2), self.work, split_interface)
        self.assertEqual(stats['attempts'], 2)
        self.assertEqual(stats['retained_triangles'], 0)
        self.assertTrue(stats['interface_recovery'])
        self.assertGreater(calls[0], 0)
        self.assertEqual(calls[-1], 0)
        validate_patch(v, t, domain, [line])
        report = audit(v, t)
        self.assertTrue(report['valid'], report)
        self.assertEqual(report['hanging_nodes'], 0)

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

    def test_raster_refinement_is_not_dropped_at_protected_constraints(self):
        # Seven-vertex reduction of stage0_r2_c0. With ordinary BAD raster
        # status, CGAL drops faces near its acute constraint cluster; the final
        # check then rejects otherwise well-shaped triangles. Use a synthetic
        # raster so the regression does not require the production DEM.
        fixture = Path(__file__).parent / 'fixtures' / 'raster_constraint_cavity.geojson'
        domain = shape(json.loads(fixture.read_text()))
        cx, cy = -1884667.712364624, 3655927.0012882375
        # Keep the production-like 50 m pixel origin fixed, including when
        # clipping the domain, since RMSE samples depend on the pixel grid.
        x0, y0 = -1885767.712364624, 3657027.0012882375
        ds = gdal.GetDriverByName('GTiff').Create(str(self.work / 'anomaly.tif'), 44, 44, 1, gdal.GDT_Float32)
        ds.SetGeoTransform((x0, 50, 0, y0, 0, -50))
        xx, yy = np.meshgrid(x0 + (np.arange(44)+.5)*50, y0 - (np.arange(44)+.5)*50)
        ds.GetRasterBand(1).WriteArray(10 + 1000*np.exp(-((xx-cx)**2 + (yy-cy)**2)/(2*300**2)))
        ds = None
        options = dict(self.options, max_area=6250000, min_area=2500,
                       dem_path=str(self.work / 'anomaly.tif'), max_tolerance=30,
                       remesh_max_vertices=5000)
        for weighted in (False, True):
            with self.subTest(weighted=weighted):
                args = dict(options, use_weights=weighted, topo_weight=1, weight_threshold=1)
                v, t = remesh_domain(domain, [], args, self.work)
                validate_patch(v, t, domain)
                self.assertTrue(audit(v, t)['valid'])
                self.assertLess(len(v), 1000)
        self.assertEqual(list(self.work.glob('cavity-*')), [])

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
                # The normal input summary is not a budget-exhaustion error.
                raise PatchError('Patch input: vertices=20, budget=10000\nPatch has quality failures')
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

    def test_refined_sloping_physical_edges_are_not_internal_interfaces(self):
        from mesher.mesh_topology import build_neighbours
        x, y = -1934381.4845238777, 3429022.4631309444
        ds = gdal.Open(str(self.dem), gdal.GA_Update)
        ds.SetGeoTransform((x-20, .5, 0, y+40, 0, -.5)); ds = None
        domain = Polygon([(x, y), (x+40.127, y+5.329), (x+40.127, y+36.329), (x, y+30.019)])
        v, t = self.mesh(domain)
        neighbours = build_neighbours(t, len(v), v)
        faces, slots = np.nonzero(neighbours < 0)
        self.assertTrue(any(not domain.boundary.covers(LineString(v[t[f, [(s+1) % 3, (s+2) % 3]], :2]))
                            for f, s in zip(faces, slots)))
        # All these triangles already belong to a conforming valid mesh. A
        # zero-size patch must not request erosion of its physical perimeter.
        joined_v, joined_t = _stitch(v, t, np.empty((0, 3)), np.empty((0, 3), dtype=int), domain)
        np.testing.assert_array_equal(joined_v[joined_t], v[t])
        validate_patch(joined_v, joined_t, domain)
        self.assertTrue(audit(joined_v, joined_t)['valid'])
        # Exercise the full merge as well, with two overlapping copies as halos.
        joined_v, joined_t, stats = merge_cavity(v, t, [x, y, x+20, y+40], v, t,
            [x+20, y, x+41, y+40], domain, [], self.options, self.work)
        self.assertGreater(stats['retained_triangles'], 0)
        validate_patch(joined_v, joined_t, domain)
        self.assertTrue(audit(joined_v, joined_t)['valid'])

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
