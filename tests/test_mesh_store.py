import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from osgeo import gdal, osr
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from mesher import mesh_store as store
from mesher.mesh_topology import build_neighbours, MeshTopologyError
from mesher.mesher_cli import build_mesh_from_arrays
from mesher.mesher_utls import MPI_do_parameterize as worker


class MeshStoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.v = np.array([[0., 0, 0], [1., 0, 0], [1., 1, 0], [0., 1, 0]])
        self.t = np.array([[0, 1, 2], [0, 2, 3]])

    def test_archive_mapping_and_disk_topology(self):
        for save in (np.savez, np.savez_compressed):
            path = self.root/'old.npz'
            save(path, verts=self.v, tris=self.t)
            with store.open_mesh_arrays(path, self.root) as (v, t):
                self.assertIsInstance(v, np.memmap)
                self.assertIsInstance(t, np.memmap)
                np.testing.assert_array_equal(v, self.v)
                neigh = build_neighbours(t, len(v), v, 1, self.root/'neigh.npy', self.root)
                self.assertIsInstance(neigh, np.memmap)
                np.testing.assert_array_equal(neigh, build_neighbours(self.t, 4))
                store.check_boundary_geometry(v, t, neigh, self.root)
            self.assertEqual(list(self.root.glob('mesh-arrays-*')), [])

    def test_disk_boundary_index_rejects_hanging_and_coincident(self):
        v = np.vstack([self.v, [[.5, 0, 0], [.5, -.5, 0], [1, -.5, 0]]])
        t = np.vstack([self.t, [[4, 5, 6]]])
        with self.assertRaisesRegex(MeshTopologyError, 'Hanging'):
            store.check_boundary_geometry(v, t, build_neighbours(t, len(v)), self.root)
        v[4] = v[0]
        with self.assertRaisesRegex(MeshTopologyError, 'Coincident'):
            store.check_boundary_geometry(v, t, build_neighbours(t, len(v)), self.root)

    def test_preparation_and_stream_exports_round_trip(self):
        ds = gdal.GetDriverByName('MEM').Create('', 2, 2, 1, gdal.GDT_Float64)
        ds.SetGeoTransform((0, 1, 0, 2, 0, -1))
        ds.GetRasterBand(1).WriteArray(np.array([[2., 3.], [4., 5.]]))
        ds.GetRasterBand(1).SetNoDataValue(-9999)
        srs = osr.SpatialReference(); srs.ImportFromEPSG(32608)
        mesh = build_mesh_from_arrays(self.v, self.t, ds, ds.GetRasterBand(1), srs, False, False, self.root)
        for key in ('vertex', 'elem', 'neigh'):
            self.assertIsInstance(mesh['mesh'][key], np.memmap)
        out = io.StringIO(); store.dump_json(mesh, out)
        decoded = json.loads(out.getvalue())
        np.testing.assert_array_equal(decoded['mesh']['elem'], self.t)
        np.testing.assert_array_equal(decoded['mesh']['vertex'], mesh['mesh']['vertex'])
        params = {'x<&': np.array([2., np.nan]), 'id': np.arange(2)}
        store.write_vtu(self.root/'test.vtu', mesh, params, {})
        reader = vtk.vtkXMLUnstructuredGridReader(); reader.SetFileName(str(self.root/'test.vtu')); reader.Update()
        result = reader.GetOutput()
        self.assertEqual(result.GetNumberOfPoints(), 4)
        self.assertEqual(result.GetNumberOfCells(), 2)
        np.testing.assert_array_equal(vtk_to_numpy(result.GetPoints().GetData()), mesh['mesh']['vertex'])
        np.testing.assert_array_equal(vtk_to_numpy(result.GetCellData().GetArray('[param] x<&')), params['x<&'])
        self.assertEqual(result.GetFieldData().GetAbstractArray('proj4').GetValue(0), srs.ExportToProj4())

    def test_parameter_worker_bounds_batches_and_writes_columns(self):
        import cloudpickle
        np.save(self.root/'v.npy', self.v)
        np.save(self.root/'t.npy', np.tile(self.t, (4, 1)))
        self.root.joinpath('config.py').write_text('mesher_debug=False\n')
        paths = {key: str(self.root/(key+'.npy')) for key in ('id', 'area')}
        for key, path in paths.items():
            store.disk_array(path, (8,), np.int64 if key == 'id' else float)
        prep = dict(nranks=1, verts_path=str(self.root/'v.npy'), elems_path=str(self.root/'t.npy'),
                    result_paths=paths, output_dir=str(self.root), batch_size=3, gt=(0,1,0,2,0,-1),
                    is_geographic=False, parameter_files={}, initial_conditions={}, RasterXSize=2,
                    RasterYSize=2, srs_proj4='+proj=utm +zone=8 +datum=WGS84', use_exactextract=True)
        inputs = self.root/'input.pickle'; inputs.write_bytes(cloudpickle.dumps(prep))
        sizes = []
        original = worker._parameterize_batch
        def batch(prep, mesh, *args):
            sizes.append(mesh['mesh']['nelem'])
            return original(prep, mesh, *args)
        with patch.object(worker, '_parameterize_batch', side_effect=batch):
            worker.main(str(inputs), False, str(self.root/'config.py'))
        self.assertEqual(sizes, [3,3,2])
        np.testing.assert_array_equal(np.load(paths['id']), np.arange(8))
        np.testing.assert_array_equal(np.load(paths['area']), np.full(8, .5))
        self.assertEqual(cloudpickle.loads((self.root/'result_0.pickle').read_bytes()),
                         dict(rank=0,start=0,stop=8,complete=True))
