import json
from itertools import product
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from shapely import covers, points, union_all
from shapely.geometry import box, LineString, mapping, Polygon, shape
from osgeo import gdal, osr

from mesher.mesher_utls import disk_merge, mpi_remesh
from mesher.mesh_topology import audit_mesh
from test_mpi_remesh import grid


@unittest.skipUnless(os.environ.get('MESHER_TEST_EXE'), 'Set MESHER_TEST_EXE')
class DiskMergeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bounds = ([0,0,40,10], [40,0,80,10])
        self.domain = box(0,0,80,10)
        self.inputs = []
        for name, bounds in zip(('a', 'b'), self.bounds):
            v, t = grid(bounds[0], bounds[2], bounds[1], bounds[3])
            path = self.root/(name+'.npz')
            np.savez(path, verts=v, tris=t)
            self.inputs.append(path)
        dem = self.root/'dem.tif'
        ds = gdal.GetDriverByName('GTiff').Create(str(dem), 82, 12, 1, gdal.GDT_Float64)
        ds.SetGeoTransform((-1, 1, 0, 11, 0, -1))
        srs = osr.SpatialReference(); srs.ImportFromEPSG(32608)
        ds.SetProjection(srs.ExportToWkt()); ds.GetRasterBand(1).Fill(10); ds = None
        self.options = dict(dem_path=str(dem), mesher_path=os.environ['MESHER_TEST_EXE'], max_area=8, min_area=.1,
                            max_tolerance=-1, errormetric='rmse', remesh_band_width=1,
                            remesh_max_attempts=6, remesh_max_vertices=10000)
        self.output = self.root/'merged.npz'

    def merge(self, **options):
        return disk_merge.merge_disk_cavity(self.inputs[0], self.bounds[0], self.inputs[1], self.bounds[1],
            self.domain, [], dict(self.options, **options), self.output)

    def assert_valid(self):
        with np.load(self.output) as mesh:
            v, t = mesh['verts'], mesh['tris']
        mpi_remesh.validate_patch(v, t, self.domain)
        neighbours = mpi_remesh.build_neighbours(t, len(v), v)
        self.assertTrue(audit_mesh({'mesh': dict(vertex=v, elem=t, neigh=neighbours,
                                                nelem=len(t), nvertex=len(v))})[0]['valid'])
        return v, t

    def test_only_seam_loaded_and_interiors_exactly_preserved(self):
        loaded = []
        original = mpi_remesh.merge_cavity
        def observe(av, at, ba, bv, bt, bb, *args, **kwargs):
            loaded.append(len(at)+len(bt))
            return original(av, at, ba, bv, bt, bb, *args, **kwargs)
        with patch.object(mpi_remesh, 'merge_cavity', side_effect=observe):
            report = self.merge()
        self.assertGreaterEqual(report['input_seconds'], 0)
        self.assertGreaterEqual(report['output_npz_seconds'], 0)
        self.assertLess(max(loaded), 400)  # 1,600 input faces; distant interiors never loaded.
        self.assertGreater(report['disk_interior_faces'], 1200)
        v, t = self.assert_valid()
        output_faces = {tuple(sorted(map(tuple, v[face]))) for face in t}
        for path in self.inputs:
            with np.load(path) as mesh:
                vv, tt = mesh['verts'], mesh['tris']
                for face in tt:
                    if np.all(vv[face, 0] < 30) or np.all(vv[face, 0] > 50):
                        self.assertIn(tuple(sorted(map(tuple, vv[face]))), output_faces)

    def test_merge_worker_publishes_numpy_face_count_as_json_integer(self):
        from mesher.mesher_utls.MPI_merge_tiles import merge_tiles
        domain = self.root / 'domain.geojson'
        domain.write_text(json.dumps(mapping(self.domain)))
        tiles = []
        for path, bounds in zip(self.inputs, self.bounds):
            marker = path.with_suffix('.json')
            marker.write_text(json.dumps({'core_bbox': bounds, 'merge_method': 'remesh'}))
            tiles.append({'npz': str(path), 'meta': str(marker)})
        count = np.count_nonzero
        with patch.object(np, 'count_nonzero', side_effect=lambda a: np.int64(count(a))):
            merge_tiles(dict(self.options, tile_a=tiles[0], tile_b=tiles[1], constraints=[],
                             outer_polygon_shp=str(domain), out_prefix=str(self.output.with_suffix(''))))
        report = json.loads(self.output.with_suffix('.json').read_text())['remesh_report']
        self.assertIs(type(report['peak_loaded_faces']), int)
        self.assertGreater(report['peak_loaded_faces'], 0)
        _, triangles = self.assert_valid()
        self.assertEqual(report['output_faces'], len(triangles))

    def test_growing_only_the_interiors_does_not_grow_the_loaded_region(self):
        small = self.merge()
        self.bounds = ([0,0,400,10], [400,0,800,10])
        self.domain = box(0,0,800,10)
        for path, bounds in zip(self.inputs, self.bounds):
            v, t = grid(bounds[0], bounds[2], bounds[1], bounds[3])
            np.savez(path, verts=v, tris=t)
        large = self.merge()
        self.assertEqual(large['peak_loaded_faces'], small['peak_loaded_faces'])
        self.assertGreater(large['disk_interior_faces'], 15000)
        self.assert_valid()

    @staticmethod
    def split_vertical_cut(v, t, fraction, x_range=(20, 40)):
        # Split one vertical boundary edge strictly inside the requested tile.
        neighbours = mpi_remesh.build_neighbours(t, len(v), v)
        for face, slot in zip(*np.nonzero(neighbours < 0)):
            a, b = t[face, [(slot+1) % 3, (slot+2) % 3]]
            if not (x_range[0] < v[a, 0] < x_range[1] and v[a, 0] == v[b, 0]):
                continue
            low, high = sorted((a, b), key=lambda i: v[i, 1])
            point = v[low] + fraction * (v[high] - v[low])
            tip = t[face, slot]
            child = len(v)
            return (np.vstack((v, point)),
                    np.vstack((np.delete(t, face, axis=0), [[tip, a, child], [tip, child, b]])), point)
        return v, t, None

    def test_split_disk_parent_matches_seam_in_one_attempt(self):
        original = mpi_remesh.merge_cavity
        points = []
        def split_boundary(*args, **kwargs):
            v, t, report = original(*args, **kwargs)
            if not points:
                # Bisect the 45 degree angle of the disk triangle.
                v, t, point = self.split_vertical_cut(v, t, np.sqrt(2) - 1)
                self.assertIsNotNone(point)
                points.append(point)
            return v, t, report

        with patch.object(mpi_remesh, 'merge_cavity', side_effect=split_boundary):
            report = self.merge()
        self.assertEqual(report['seam_attempts'], 1)
        self.assertEqual(report['disk_parents_subdivided'], 1)
        self.assertEqual(report['unmatched_disk_edges'], 1)
        self.assertGreater(report['disk_interior_faces'], 1200)
        v, _ = self.assert_valid()
        self.assertEqual(np.count_nonzero((v == points[0]).all(axis=1)), 1)

    def test_changed_outer_interface_expands_and_reuses_generated_work(self):
        original = mpi_remesh.merge_cavity
        calls = []
        def split_boundary(*args, **kwargs):
            calls.append(kwargs.get('initial_keep'))
            v, t, report = original(*args, **kwargs)
            if len(calls) == 1:
                # The midpoint would create a disk child below the angle bound.
                v, t, point = self.split_vertical_cut(v, t, .5)
                self.assertIsNotNone(point)
            return v, t, report
        with patch.object(mpi_remesh, 'merge_cavity', side_effect=split_boundary):
            report = self.merge()
        self.assertEqual(len(calls), 2)
        self.assertEqual(report['attempts'], 1)
        self.assertEqual(report['original_triangles_retained'], np.count_nonzero(calls[1]))
        self.assertGreater(report['seam_work_reused'], 0)
        self.assertGreater(report['disk_interior_faces'], 1000)
        self.assert_valid()

    def test_mixed_disk_subdivisions_expand_only_the_failed_parent(self):
        # Both interiors exceed the spare face budget. One exact split fits its
        # disk parent and one would create an 18.4-degree child. Only the
        # failed parent needs to draw its neighbours into the next solve.
        original = mpi_remesh.merge_cavity
        splits = []
        def split_both_cuts(av, at, *args, **kwargs):
            if splits:
                kept_xy = {tuple(xy) for xy in av[np.unique(at[kwargs['initial_keep']]), :2]}
                self.assertIn(tuple(splits[0][:2]), kept_xy)
                # The right-hand parent has been repaired, so its disk-side
                # neighbour collar must still be on disk.
                right_cut = splits[0][0]
                self.assertLessEqual(av[:, 0].max(), right_cut + 1)
            v, t, report = original(av, at, *args, **kwargs)
            if not splits:
                v, t, good = self.split_vertical_cut(v, t, 2-np.sqrt(2), (40, 60))
                v, t, bad = self.split_vertical_cut(v, t, .5, (20, 40))
                self.assertIsNotNone(good)
                self.assertIsNotNone(bad)
                splits.extend((good, bad))
            return v, t, report

        with patch.object(mpi_remesh, 'merge_cavity', side_effect=split_both_cuts):
            report = self.merge(remesh_max_attempts=2, remesh_max_active_faces=500)
        self.assertEqual(report['seam_attempts'], 2)
        first, = report['seam_attempt_reports']
        self.assertEqual(first['unmatched_disk_edges'], 2)
        self.assertEqual(first['disk_parents_subdivided'], 1)
        self.assertEqual(first['unresolved_disk_edges'], 1)
        self.assertEqual(report['unresolved_disk_edges'], 0)
        self.assertGreater(report['disk_interior_faces'], 1000)
        self.assertLessEqual(report['peak_loaded_faces'], 500)
        v, _ = self.assert_valid()
        self.assertEqual(np.count_nonzero((v == splits[0]).all(axis=1)), 1)

    def moving_cut(self, splits):
        original = mpi_remesh.merge_cavity
        def split_boundary(*args, **kwargs):
            v, t, report = original(*args, **kwargs)
            v, t, point = self.split_vertical_cut(v, t, .5)
            splits.append(point is not None)
            return v, t, report
        return patch.object(mpi_remesh, 'merge_cavity', side_effect=split_boundary)

    def test_stalled_cut_clears_refined_band_without_loading_whole_interior(self):
        # Model a backend that splits any exposed vertical edge in a fixed
        # strip. The prescribed grid elsewhere keeps this moving-cut regression
        # independent of CGAL refinement order. Disk subdivisions use CGAL.
        vertices, triangles = grid(0, 80, 0, 10)
        centers = points(vertices[triangles, :2].mean(axis=1))
        splits = []
        def split_inside_strip(*args, prepared_domain, **kwargs):
            faces = triangles[covers(prepared_domain, centers)]
            used, index = np.unique(faces, return_inverse=True)
            v, t, point = self.split_vertical_cut(vertices[used], index.reshape(-1, 3), .5, (32, 40))
            splits.append(point is not None)
            return v, t, {}

        with patch.object(mpi_remesh, 'merge_cavity', side_effect=split_inside_strip):
            report = self.merge(remesh_max_attempts=4, remesh_max_active_faces=500)
        self.assertEqual(splits, [True, True, True, False])
        self.assertLessEqual(report['peak_loaded_faces'], 500)
        self.assertGreater(report['disk_interior_faces'], 1000)
        self.assert_valid()

    def test_repeated_moving_cut_closes_only_affected_interior_within_budget(self):
        splits = []
        with self.moving_cut(splits):
            report = self.merge(remesh_max_attempts=3, remesh_max_active_faces=1000)
        self.assertEqual(splits, [True, True, False])
        self.assertEqual(report['seam_attempts'], 3)
        self.assertTrue(report['seam_attempt_reports'][-1]['closing_disk_components'])
        self.assertLessEqual(report['peak_loaded_faces'], 1000)
        self.assertGreater(report['disk_interior_faces'], 600)
        self.assert_valid()

    def test_final_cached_repair_handles_a_new_cut_on_the_other_tile(self):
        # A cached candidate split the left disk cut. On its last solve CGAL
        # can instead split the previously matching right cut. Both collars
        # fit, but an unrelated island on the right must stay on disk.
        gap = box(60, 0, 61, 10)
        self.domain = self.domain.difference(gap)
        with np.load(self.inputs[1]) as tile:
            v, t = tile['verts'], tile['tris']
        x = v[t, 0].mean(axis=1)
        np.savez(self.inputs[1], verts=v, tris=t[(x < 60) | (x > 61)])
        island_faces = int((x > 61).sum())

        original = mpi_remesh.merge_cavity
        splits = []
        def split_opposite_sides(*args, **kwargs):
            v, t, report = original(*args, **kwargs)
            side = (20, 40) if not splits else (40, 60)
            v, t, point = self.split_vertical_cut(v, t, .5, side)
            splits.append(point is not None)
            return v, t, report

        with patch.object(mpi_remesh, 'merge_cavity', side_effect=split_opposite_sides):
            report = self.merge(remesh_max_attempts=2, remesh_max_active_faces=1300)
        self.assertEqual(splits, [True, False])
        self.assertEqual(report['seam_attempts'], 2)
        self.assertGreater(report['seam_work_reused'], 0)
        self.assertLessEqual(report['peak_loaded_faces'], 1300)
        self.assertEqual(report['disk_interior_faces'], island_faces)
        self.assert_valid()

    def test_disconnected_interior_on_same_tile_does_not_consume_recovery_budget(self):
        gap = box(19, 0, 20, 10)
        self.domain = self.domain.difference(gap)
        with np.load(self.inputs[0]) as tile:
            v, t = tile['verts'], tile['tris']
        x = v[t, 0].mean(axis=1)
        np.savez(self.inputs[0], verts=v, tris=t[(x < 19) | (x > 20)])

        with self.moving_cut([]):
            report = self.merge(remesh_max_attempts=3, remesh_max_active_faces=530)
        self.assertTrue(report['seam_attempt_reports'][-1]['closing_disk_components'])
        self.assertLessEqual(report['peak_loaded_faces'], 530)
        self.assertGreater(report['disk_interior_faces'], 1100)
        self.assert_valid()

    def test_final_solve_after_rejected_candidate_closes_disk_cuts(self):
        # A rejection leaves no cached seam to extend. The final solve must
        # still close affordable cuts even if the new local mesh splits them.
        original = mpi_remesh.merge_cavity
        coverage = disk_merge.inactive_coverage
        cuts, loaded = {}, []
        def observe_cut(vertices, *args):
            result = coverage(vertices, *args)
            cuts[id(vertices)] = result[1]
            return result

        def reject_then_split(av, at, ba, bv, bt, bb, domain, constraints, options, directory, **kwargs):
            loaded.append(len(at)+len(bt))
            if len(loaded) == 1:
                raise mpi_remesh.PatchError('Local cavity rejected')
            for edges in cuts.values():
                if edges:
                    # Require a short segment of an actual disk edge. The
                    # native solver produces a valid local mesh, but the disk
                    # parent cannot absorb this split at the required angle.
                    a, b = min(edges)
                    near = np.asarray(a) + (np.asarray(b)-a)/32
                    segment = LineString([a, near, b])
                    region = kwargs['prepared_domain']
                    v, t = mpi_remesh.remesh_domain(region, [segment], options, directory)
                    return v, t, {}
            return original(av, at, ba, bv, bt, bb, domain, constraints, options, directory, **kwargs)

        with patch.object(mpi_remesh, 'merge_cavity', side_effect=reject_then_split), \
             patch.object(disk_merge, 'inactive_coverage', side_effect=observe_cut):
            report = self.merge(remesh_max_attempts=2, remesh_max_active_faces=1600)
        self.assertEqual(loaded[-1], 1600)
        self.assertLess(loaded[0], 400)
        self.assertEqual(report['seam_attempts'], 2)
        self.assertEqual(report['disk_interior_faces'], 0)
        self.assert_valid()

    def test_exhausted_moving_cut_names_inputs_and_preserves_previous_output(self):
        splits = []
        self.output.write_bytes(b'previous completed mesh')
        with self.moving_cut(splits), self.assertRaises(mpi_remesh.PatchError) as failure:
            self.merge(remesh_max_attempts=3, remesh_max_active_faces=500)
        self.assertRegex(str(failure.exception), r'Seam-region expansion exhausted.*merged.*unmatched disk edges')
        for path in self.inputs:
            self.assertIn(str(path), str(failure.exception))
        self.assertIn('active_face_limit=500', str(failure.exception))
        self.assertEqual(splits, [True, True, True])
        self.assertEqual(self.output.read_bytes(), b'previous completed mesh')

    def test_active_budget_is_checked_before_materialization_and_keeps_previous_output(self):
        self.output.write_bytes(b'previous completed mesh')
        with patch.object(disk_merge, 'compact_region', side_effect=AssertionError('loaded seam')):
            with self.assertRaisesRegex(mpi_remesh.PatchError, 'active-face budget'):
                self.merge(remesh_max_active_faces=1)
        self.assertEqual(self.output.read_bytes(), b'previous completed mesh')
        self.assertFalse(list(self.root.glob('seam-store-*')))

    def test_native_failure_is_not_retried_or_published(self):
        self.output.write_bytes(b'previous completed mesh')
        with patch.object(mpi_remesh, 'merge_cavity', side_effect=mpi_remesh.BackendTerminated('signal 11')) as native:
            with self.assertRaises(mpi_remesh.BackendTerminated):
                self.merge()
        self.assertEqual(native.call_count, 1)
        self.assertEqual(self.output.read_bytes(), b'previous completed mesh')

    def test_empty_pair_remains_empty(self):
        for path in self.inputs:
            np.savez(path, verts=np.empty((0,3)), tris=np.empty((0,3), dtype=np.int64))
        self.domain = box(100,100,110,110)
        report = self.merge()
        self.assertEqual(report['output_faces'], 0)
        with np.load(self.output) as output:
            self.assertEqual(output['verts'].shape, (0,3))
            self.assertEqual(output['tris'].shape, (0,3))

    def test_final_topology_failure_cannot_replace_existing_archive(self):
        from mesher.mesh_topology import MeshTopologyError
        original = disk_merge.stream_join
        def corrupt(*args, **kwargs):
            v, t = original(*args, **kwargs)
            t[-1] = t[0]
            return v, t
        self.output.write_bytes(b'previous completed mesh')
        with patch.object(disk_merge, 'stream_join', side_effect=corrupt):
            with self.assertRaises(MeshTopologyError):
                self.merge()
        self.assertEqual(self.output.read_bytes(), b'previous completed mesh')

    def test_production_roundoff_components_release_only_adjacent_faces(self):
        path = Path(__file__).parent / 'fixtures' / 'roundoff_boundary_cavity.geojson'
        cavity = shape(json.loads(path.read_text()))
        parts = list(cavity.geoms)
        self.check_roundoff_merge(cavity, parts[:3], parts[3:], box(-1917200, 3717000, -1917000, 3717200),
            ([-1950000,3680000,-1890000,3740000], [-1890000,3680000,-1880000,3740000]), 20)

    def test_coarse_tvc_roundoff_before_only_permitted_solve(self):
        path = Path(__file__).parent / 'fixtures' / 'tvc_coarse_roundoff.geojson'
        cavity = shape(json.loads(path.read_text()))
        main, *slivers = cavity.geoms
        self.assertAlmostEqual(cavity.area, 992000646.605899, places=5)
        bounds = ([-1524150.8290839698, 3516737.591255273, -1420000, 3574638.1809423035],
                  [-1420000, 3516737.591255273, -1400819.7600370413, 3574638.1809423035])
        for tolerance in (0.25, 30, 500):
            with self.subTest(tolerance_m=tolerance):
                self.check_roundoff_merge(cavity, slivers, [main], box(-1510000, 3560000, -1509800, 3560200),
                                          bounds, tolerance)

    def check_roundoff_merge(self, cavity, slivers, resolved_parts, island, bounds, tolerance):
        self.assertEqual(len(mpi_remesh._roundoff_chords(cavity)), len(slivers))

        # Reconstruct the unavailable source triangles on one side of each
        # supplied chord. Subtracting these fans reproduces the exact cavity.
        vertices, triangles, parents = [], [], []
        for part in slivers:
            a, b, *chain = np.asarray(part.exterior.coords[:-1])
            normal = np.array([a[1]-b[1], b[0]-a[0]])
            apex = (a+b)/2 + np.sign(normal @ (chain[0]-a)) * .3 * normal
            parents.append(Polygon([a, b, apex]))
            start = len(vertices)
            vertices.extend([a, *reversed(chain), b, apex])
            tip = len(vertices)-1
            triangles.extend([i, i+1, tip] for i in range(start, tip-1))

        nparents = len(triangles)
        start = len(vertices)
        vertices.extend(island.exterior.coords[:-1])
        triangles.extend([[start, start+1, start+2], [start, start+2, start+3]])
        xy, ids = np.unique(vertices, axis=0, return_inverse=True)
        v = np.column_stack((xy, np.zeros(len(xy))))
        t = ids[np.asarray(triangles)]
        self.domain = union_all([*parents, *resolved_parts, island])
        self.assertTrue(self.domain.difference(mpi_remesh.mesh_coverage(v, t)).equals(cavity))
        self.bounds = bounds
        empty_v, empty_t = v[:0], t[:0]
        np.savez(self.inputs[0], verts=v, tris=t)
        np.savez(self.inputs[1], verts=empty_v, tris=empty_t)
        x0, y0, x1, y1 = self.domain.bounds
        ds = gdal.GetDriverByName('GTiff').Create(str(self.root/'flat.tif'),
            int((x1-x0)/100)+3, int((y1-y0)/100)+3, 1, gdal.GDT_Float64)
        ds.SetGeoTransform((x0-100, 100, 0, y1+100, 0, -100))
        ds.GetRasterBand(1).Fill(10)
        ds = None
        options = dict(self.options, dem_path=str(self.root/'flat.tif'), max_area=6250000,
                       min_area=2500, max_tolerance=tolerance, remesh_band_width=1000, remesh_max_vertices=500000,
                       remesh_max_attempts=1)

        # These faces begin on disk. Planning must enforce the budget before
        # loading any released parent, even with only one solve allowed.
        self.output.write_bytes(b'previous completed mesh')
        with patch.object(disk_merge, 'compact_region', side_effect=AssertionError('loaded seam')):
            with self.assertRaisesRegex(mpi_remesh.PatchError, 'active-face budget'):
                self.merge(**options, remesh_max_active_faces=nparents-1)
        self.assertEqual(self.output.read_bytes(), b'previous completed mesh')

        for disk in (False, True):
            with self.subTest(disk=disk), patch.object(mpi_remesh, '_run_backend',
                                                       wraps=mpi_remesh._run_backend) as backend:
                if disk:
                    report = self.merge(**options)
                    result_v, result_t = self.assert_valid()
                    self.assertEqual(report['disk_interior_faces'], 2)
                    self.assertEqual(report['roundoff_disk_parents'], nparents)
                    self.assertEqual(report['attempts'], 1)
                    self.assertEqual(report['seam_attempts'], 1)
                else:
                    result_v, result_t, report = mpi_remesh.merge_cavity(
                        v, t, self.bounds[0], empty_v, empty_t, self.bounds[1],
                        self.domain, [], options, self.root)
                    self.assertEqual(report['original_triangles_retained'], 2)
                    self.assertEqual(report['attempt_stats'][0]['roundoff_parents_released'], nparents)
                    self.assertEqual(report['attempts'], 1)
                self.assertEqual(backend.call_count, 1)
                mpi_remesh.validate_patch(result_v, result_t, self.domain)
                output = {tuple(sorted(map(tuple, result_v[face]))) for face in result_t}
                for face in t[-2:]:
                    self.assertIn(tuple(sorted(map(tuple, v[face]))), output)

    def test_unbisectable_corner_releases_its_parent_before_native_refinement(self):
        cavity = shape(json.loads((Path(__file__).parent/'fixtures'/'tvc_r2c1_unbisectable_corner.geojson').read_text()))
        a, b, c, *others = cavity.exterior.coords[:-1]
        self.assertEqual(np.nextafter(a[0], b[0]), b[0])
        self.assertEqual(a[1], b[1])
        # Reconstruct a parent whose subtraction leaves the exact reported
        # 11-vertex cavity, including its one-coordinate-increment bevel.
        far = (b[0]+3000, b[1])
        self.domain = Polygon([a, far, c, *others])
        parent = Polygon([b, far, c])
        self.assertTrue(self.domain.difference(parent).equals(cavity))
        v = np.column_stack((parent.exterior.coords[:-1], np.zeros(3)))
        t = np.array([[0, 1, 2]])
        self.bounds = ([-1821000, 3667000, -1810000, 3678000], [-1810000, 3667000, -1800000, 3678000])
        np.savez(self.inputs[0], verts=v, tris=t)
        np.savez(self.inputs[1], verts=v[:0], tris=t[:0])
        ds = gdal.Open(self.options['dem_path'], gdal.GA_Update)
        ds.SetGeoTransform((-1821000, 1000, 0, 3678000, 0, -1000))
        ds = None
        for tolerance, disk in product((0.25, 20, 500), (False, True)):
            opts = dict(self.options, max_area=6250000, min_area=2500, max_tolerance=tolerance,
                        remesh_max_attempts=1)
            with self.subTest(tolerance_m=tolerance, disk=disk), patch.object(
                    mpi_remesh, '_run_backend', wraps=mpi_remesh._run_backend) as backend:
                if disk:
                    report = self.merge(**opts)
                    self.assertEqual(report['roundoff_disk_parents'], 1)
                    result_v, result_t = self.assert_valid()
                else:
                    result_v, result_t, report = mpi_remesh.merge_cavity(v, t, self.bounds[0],
                        v[:0], t[:0], self.bounds[1], self.domain, [], opts, self.root)
                    self.assertEqual(report['attempt_stats'][0]['roundoff_parents_released'], 1)
                self.assertEqual(backend.call_count, 1)
                self.assertEqual(report['attempts'], 1)
                mpi_remesh.validate_patch(result_v, result_t, self.domain)

    def test_physical_hole_in_disk_interior_is_preserved(self):
        hole = box(5, 3, 7, 5)
        self.domain = self.domain.difference(hole)
        path = self.inputs[0]
        with np.load(path) as mesh:
            v, t = mesh['verts'], mesh['tris']
        centers = v[t, :2].mean(axis=1)
        t = t[~((centers[:,0]>5) & (centers[:,0]<7) & (centers[:,1]>3) & (centers[:,1]<5))]
        np.savez(path, verts=v, tris=t)
        report = self.merge()
        self.assertGreater(report['disk_interior_faces'], 1200)
        self.assert_valid()
