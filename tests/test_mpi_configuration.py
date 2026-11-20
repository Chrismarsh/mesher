"""Migration and reuse rules for cavity-only MPI meshing."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from mesher.mesher_cli import (mpi_remesh_options_from_config, validate_mpi_mesh_metadata,
                               compute_tile_grid, compute_mpi_final_npz_path, mpi_final_npz_name)
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

    def test_disabled_smoothing_setting_is_compatible(self):
        defaults = mpi_remesh_options_from_config(SimpleNamespace())
        self.assertEqual(defaults, mpi_remesh_options_from_config(SimpleNamespace(mpi_global_lloyd=0)))

    def test_retired_smoother_controls_fail_clearly(self):
        for count in (1, -1):
            with self.subTest(count=count), self.assertRaisesRegex(ValueError, 'smoothing.*retired.*mpi_global_lloyd'):
                mpi_remesh_options_from_config(SimpleNamespace(mpi_global_lloyd=count))
        with self.assertRaisesRegex(ValueError, 'mpi_lloyd_boundary_tol.*retired.*remove'):
            mpi_remesh_options_from_config(SimpleNamespace(mpi_lloyd_boundary_tol=0))

    def test_mpi_cannot_enable_serial_lloyd(self):
        with self.assertRaisesRegex(ValueError, 'lloyd_itr must be 0 with MPI'):
            mpi_remesh_options_from_config(SimpleNamespace(lloyd_itr=1))


class MpiTilingTests(unittest.TestCase):
    def test_allocations_no_longer_round_down_to_square(self):
        self.assertEqual(compute_tile_grid(95, 1, 1), (9, 10))
        self.assertEqual(compute_tile_grid(96, 1, 1), (8, 12))
        self.assertEqual(compute_tile_grid(96, 616655, 521105), (8, 12))

    def test_physical_dimensions_determine_orientation(self):
        self.assertEqual(compute_tile_grid(16, 4, 1), (2, 8))
        self.assertEqual(compute_tile_grid(16, 1, 4), (8, 2))
        for count in (1, 2, 3, 5, 7, 11, 50, 95, 96, 97):
            with self.subTest(count=count):
                rows, cols = compute_tile_grid(count, 4, 1)
                self.assertGreaterEqual(rows * cols, 1)
                self.assertLessEqual(rows * cols, count)
        self.assertEqual(compute_tile_grid(7, 1, 1), (2, 3))

    def test_invalid_grid_inputs_fail_before_mpi_launch(self):
        for args in ((0, 1, 1), (-1, 1, 1), (True, 1, 1), (2.5, 1, 1),
                     (4, 0, 1), (4, 1, float('nan')), (4, float('inf'), 1)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                compute_tile_grid(*args)

    def test_manifest_reuse_does_not_recompute_layout_from_new_rank_count(self):
        with tempfile.TemporaryDirectory() as directory:
            final = mpi_final_npz_name('case', 8, 12)
            path = Path(directory) / 'case_mpi_layout.json'
            path.write_text(json.dumps({'version': 1, 'rows': 8, 'cols': 12,
                                        'mpi_nworkers': 96, 'final_npz': final}))
            self.assertEqual(compute_mpi_final_npz_path(directory, 'case', 50),
                             str(Path(directory) / final))

    def test_legacy_completed_square_layout_still_resolves(self):
        with tempfile.TemporaryDirectory() as directory:
            for ranks, stage in ((50, 5), (95, 7), (96, 7)):
                with self.subTest(ranks=ranks):
                    self.assertEqual(compute_mpi_final_npz_path(directory, 'case', ranks),
                                     str(Path(directory) / f'case_mpi_stage{stage}_r0_c0.npz'))
            self.assertEqual(compute_mpi_final_npz_path(directory, 'case', 1),
                             str(Path(directory) / 'case_mpi_tile_0.npz'))

    def test_invalid_manifest_is_not_silently_replaced_with_legacy_guess(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'case_mpi_layout.json'
            for content in ('not json', '[]', '{}',
                            json.dumps({'version': 1, 'rows': 2, 'cols': 3, 'final_npz': '../wrong.npz'})):
                path.write_text(content)
                with self.subTest(content=content), self.assertRaisesRegex(ValueError, 'regenerate'):
                    compute_mpi_final_npz_path(directory, 'case', 96)
