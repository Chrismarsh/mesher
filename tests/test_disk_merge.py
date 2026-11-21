import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from shapely.geometry import box, mapping
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

    def test_changed_outer_interface_expands_and_reuses_generated_work(self):
        original = mpi_remesh.merge_cavity
        calls = []
        def split_boundary(*args, **kwargs):
            calls.append(kwargs.get('initial_keep'))
            v, t, report = original(*args, **kwargs)
            if len(calls) == 1:
                neighbours = mpi_remesh.build_neighbours(t, len(v), v)
                faces, slots = np.nonzero(neighbours < 0)
                for f, s in zip(faces, slots):
                    a, b = t[f, [(s+1)%3, (s+2)%3]]
                    if np.all((v[[a,b], 0] > 1) & (v[[a,b], 0] < 79)) and v[a,0] == v[b,0]:
                        midpoint = (v[a]+v[b])/2
                        nv = len(v)
                        v = np.vstack((v, midpoint))
                        t = np.vstack((np.delete(t, f, axis=0), [[t[f,s], a, nv], [t[f,s], nv, b]]))
                        break
                else:
                    self.fail('fixture has no disk interface')
            return v, t, report
        with patch.object(mpi_remesh, 'merge_cavity', side_effect=split_boundary):
            report = self.merge()
        self.assertEqual(len(calls), 2)
        self.assertEqual(report['attempts'], 1)
        self.assertEqual(report['original_triangles_retained'], np.count_nonzero(calls[1]))
        self.assertGreater(report['seam_work_reused'], 0)
        self.assertGreater(report['disk_interior_faces'], 1000)
        self.assert_valid()

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
