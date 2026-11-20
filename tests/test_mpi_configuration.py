"""Migration and reuse rules for cavity-only MPI meshing."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from mesher.mesher_cli import mpi_remesh_options_from_config, validate_mpi_mesh_metadata
from mesher.mesher_utls.MPI_mesh_tiles import mesh_tile
from mesher.mesher_utls.MPI_merge_tiles import merge_tiles


class MpiConfigurationTests(unittest.TestCase):
    def test_remesh_is_default_and_explicit_choice(self):
        defaults = mpi_remesh_options_from_config(SimpleNamespace())
        self.assertEqual(defaults['merge_method'], 'remesh')
        self.assertEqual(defaults, mpi_remesh_options_from_config(SimpleNamespace(mpi_merge_method='remesh')))
        configured = mpi_remesh_options_from_config(SimpleNamespace(
            mpi_remesh_band_width=100, mpi_remesh_max_attempts=5, mpi_remesh_max_vertices=2000))
        self.assertEqual(configured, {'merge_method': 'remesh', 'remesh_band_width': 100,
                                     'remesh_max_attempts': 5, 'remesh_max_vertices': 2000})

    def test_retired_method_and_controls_fail_with_migration_message(self):
        for method in ('shared', 'typo'):
            with self.subTest(method=method), self.assertRaisesRegex(ValueError, 'Only.*remesh.*retired'):
                mpi_remesh_options_from_config(SimpleNamespace(mpi_merge_method=method))
        for key, value in (('mpi_merge_snap_tol', 1e-6), ('mpi_shared_edge_spacing', 100),
                           ('mpi_shared_edge_constraints', True)):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, key + '.*remove it'):
                mpi_remesh_options_from_config(SimpleNamespace(**{key: value}))

    def test_tile_and_merge_workers_reject_old_tasks_before_reading_inputs(self):
        for worker in (mesh_tile, merge_tiles):
            with self.subTest(worker=worker.__name__), self.assertRaisesRegex(ValueError, 'Shared-edge.*retired'):
                worker({'merge_method': 'shared'})

    def test_completed_remesh_metadata_can_be_reused_without_modification(self):
        with tempfile.TemporaryDirectory() as directory:
            mesh = Path(directory) / 'joined.npz'
            meta = mesh.with_suffix('.json')
            mesh.write_bytes(b'unchanged mesh')
            contents = json.dumps({'merge_method': 'remesh', 'core_bbox': [0, 0, 10, 10]})
            meta.write_text(contents)
            validate_mpi_mesh_metadata(mesh)
            self.assertEqual(mesh.read_bytes(), b'unchanged mesh')
            self.assertEqual(meta.read_text(), contents)

    def test_shared_missing_and_invalid_metadata_cannot_be_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            mesh = Path(directory) / 'joined.npz'
            meta = mesh.with_suffix('.json')
            mesh.write_bytes(b'unchanged mesh')
            for contents in (None, '{}', '{"merge_method":"shared"}', 'invalid json', '[]'):
                if contents is not None:
                    meta.write_text(contents)
                with self.subTest(contents=contents), self.assertRaisesRegex(ValueError, 'reuse_mesh=False'):
                    validate_mpi_mesh_metadata(mesh)
                self.assertEqual(mesh.read_bytes(), b'unchanged mesh')
                if contents is not None:
                    self.assertEqual(meta.read_text(), contents)
                    meta.unlink()

    def test_merge_worker_refuses_shared_tile_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            meta = Path(directory) / 'old.json'
            meta.write_text(json.dumps({'core_bbox': [0, 0, 10, 10]}))
            tile = {'meta': str(meta), 'npz': 'not-read.npz'}
            with self.assertRaisesRegex(ValueError, 'Cannot merge old shared-edge tiles'):
                merge_tiles({'tile_a': tile, 'tile_b': tile})
