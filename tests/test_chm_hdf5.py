import importlib.util
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from mesher.chm_hdf5 import export_hdf5
from mesher.mesh_store import batches, copy_array
from mesher.mesh_topology import build_neighbours


@unittest.skipUnless(importlib.util.find_spec('h5py'), 'h5py is optional')
class ChmHdf5Tests(unittest.TestCase):
    # Use file-backed geometry and recognizable per-face values so a misplaced row can be traced
    # through the export.
    def setUp(self):
        import h5py
        self.h5py = h5py
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.directory = Path(tmp.name)
        self.prefix = self.directory / 'terrain'
        self.vertices = np.array([[x, y, x+y] for y in range(3) for x in range(3)], dtype=float)
        self.faces = np.array([[a, a+1, a+4] if half == 0 else [a, a+4, a+3]
                               for a in (0, 1, 3, 4) for half in range(2)])
        self.neighbours = build_neighbours(self.faces, 9)
        self.mesh = {k: copy_array(self.directory / (k+'.npy'), v) for k, v in
                     [('vertex', self.vertices), ('elem', self.faces), ('neigh', self.neighbours)]}
        self.mesh['is_geographic'] = True
        self.params = {k: copy_array(self.directory / (k+'.npy'), v) for k, v in
                       [('id', np.arange(8)), ('snow', np.arange(8, dtype=float) + .125)]}
        self.proj4 = '+proj=longlat +datum=WGS84 +no_defs'

    # Read the files the way CHM expects, checking array datatypes, root attributes, neighbour
    # numbering and parameter correspondence under the supplied order.
    def check_files(self, sizes, order):
        with self.h5py.File(self.directory / 'terrain_mesh.h5') as mesh:
            g = mesh['mesh']
            for name, count, dt in [('vertex', 9, 'f8'), ('elem', 8, 'i4'), ('neighbor', 8, 'i4')]:
                self.assertEqual(g[name].shape, (count,))
                self.assertEqual(g[name].dtype.subdtype, (np.dtype(dt), (3,)))
            np.testing.assert_equal(g['vertex'][:], self.vertices)
            np.testing.assert_equal(g['elem'][:], self.faces[order])
            inverse = np.argsort(order)
            old = self.neighbours[order]
            np.testing.assert_equal(g['neighbor'][:], np.where(old == -1, -1, inverse[old]))
            np.testing.assert_equal(g['local_sizes'][:], sizes)
            np.testing.assert_equal(g['cell_global_id'][:], np.arange(8))
            self.assertEqual(mesh.attrs['/mesh/proj4'][0], self.proj4.encode())
            self.assertEqual(mesh.attrs['/mesh/partition_method'][0], b'metis' if len(sizes) > 1 else b'')
            self.assertEqual(mesh.attrs['/mesh/version'][0], b'2.0.0' if len(sizes) > 1 else b'1.0.0')
            self.assertEqual(mesh.attrs['/mesh/is_partition'][0], 0)
            self.assertEqual(mesh.attrs['/mesh/is_geographic'][0], 1)
        with self.h5py.File(self.directory / 'terrain_param.h5') as params:
            for key, values in self.params.items():
                dataset = params['parameters'][key]
                self.assertEqual(dataset.dtype, np.float64)
                np.testing.assert_equal(dataset[:], values[order])

    # Force short writes across permutation boundaries so an incorrect inverse neighbour map
    # cannot hide within a single batch.
    def test_serial_identity_and_nontrivial_order(self):
        export_hdf5(self.prefix, self.mesh, self.params, self.proj4, 1, 1, None)
        self.check_files([8], np.arange(8))
        order = np.array([7, 4, 1, 5, 2, 0, 3, 6])
        with patch('mesher.tools.meshpermutation.metis_face_order', return_value=(np.array([3, 5]), order)), \
                patch('mesher.chm_hdf5.batches', side_effect=lambda count: batches(count, 2)):
            export_hdf5(self.prefix, self.mesh, self.params, self.proj4, 2, 1, None)
        self.check_files([3, 5], order)

    # Simulate running out of disk space during a write and check that a reader still sees the
    # previous complete pair.
    def test_failed_export_preserves_prior_outputs(self):
        export_hdf5(self.prefix, self.mesh, self.params, self.proj4, 1, 1, None)
        with patch('mesher.chm_hdf5._write_slices', side_effect=OSError('disk full')):
            with self.assertRaisesRegex(OSError, 'disk full'):
                export_hdf5(self.prefix, self.mesh, self.params, self.proj4, 1, 1, None)
        self.check_files([8], np.arange(8))

    # Read back a file produced by real MPI writers to catch collective schema or slice-order
    # errors that a serial export cannot exercise.
    @unittest.skipUnless(shutil.which('mpirun') and os.environ.get('MESHER_TEST_MPI_HDF5'),
                         'enable with MESHER_TEST_MPI_HDF5=1 on a parallel h5py installation')
    def test_parallel_metis_export(self):
        sizes = np.array([3, 5])
        order = np.array([7, 4, 1, 5, 2, 0, 3, 6])
        with patch('mesher.tools.meshpermutation.metis_face_order', return_value=(sizes, order)):
            export_hdf5(self.prefix, self.mesh, self.params, self.proj4, 2, 3, f'mpirun -n 3 {sys.executable}')
        self.check_files(sizes, order)

    # Reuse a completed mesh to exercise the CLI's HDF5 and optional text outputs independently
    # of terrain refinement.
    @unittest.skipUnless(importlib.util.find_spec('exactextract') and shutil.which('mpirun')
                         and shutil.which('gdal-config'), 'requires the Mesher CLI dependencies')
    def test_cli_export_switches_with_reused_mpi_mesh(self):
        from osgeo import gdal, osr
        import subprocess
        env = os.environ.copy()
        env['PYTHONPATH'] = os.pathsep.join(str(Path(entry).resolve()) for entry in sys.path)
        dem = self.directory / 'dem.tif'
        crs = osr.SpatialReference()
        crs.ImportFromEPSG(32608)
        ds = gdal.GetDriverByName('GTiff').Create(str(dem), 6, 6, 1, gdal.GDT_Float32)
        ds.SetGeoTransform((500000, 1, 0, 6500006, 0, -1))
        ds.SetProjection(crs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(np.ones((6, 6)))
        ds.GetRasterBand(1).SetNoDataValue(-9999)
        ds = None
        for mode, settings in [('legacy', ''), ('binary', 'write_hdf5=True\nwrite_text=False\nwrite_param_per_file=True'),
                               ('both', 'write_hdf5=True\nwrite_param_per_file=True')]:
            with self.subTest(mode=mode):
                out = self.directory / mode
                work = out / 'dem'
                work.mkdir(parents=True)
                np.savez(work / 'dem_mpi_tile_0.npz',
                         verts=self.vertices[:, :2] + [500001, 6500001], tris=self.faces)
                (work / 'dem_mpi_tile_0.json').write_text('{"merge_method":"remesh"}')
                cfg = self.directory / (mode + '.py')
                cfg.write_text(f"dem_filename='dem.tif'\nmax_area=10\nmpi_mesh=True\nreuse_mesh=True\n"
                               f"user_output_dir='{mode}'\nwrite_shp=False\nwrite_vtu=False\nMPI_nworkers=1\n"
                               f"MPI_exec_str='mpirun -n 1 {sys.executable}'\n{settings}\n")
                result = subprocess.run([sys.executable, '-m', 'mesher.mesher_cli', str(cfg)],
                                        cwd=self.directory, env=env, capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual((out / 'dem_mpi.mesh').exists(), mode != 'binary')
                self.assertEqual((out / 'dem_mpi.param').exists(), mode == 'legacy')
                self.assertEqual(sorted(p.name for p in out.glob('*.param.*.json')),
                                 ['dem_mpi.param.area.json', 'dem_mpi.param.id.json'] if mode == 'both' else [])
                self.assertEqual((out / 'dem_mpi_mesh.h5').exists(), mode != 'legacy')
                self.assertEqual((out / 'dem_mpi_param.h5').exists(), mode != 'legacy')
                self.assertTrue((out / 'dem_mpi.ic').exists())
