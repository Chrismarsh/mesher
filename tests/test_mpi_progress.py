import io
import json
from pathlib import Path
import shlex
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch, MagicMock
from types import SimpleNamespace
import subprocess

import cloudpickle

from mesher.mesher_utls.mpi_progress import StageProgress, complete_task


class Terminal(io.StringIO):
    def isatty(self):
        return True


class MpiProgressTests(unittest.TestCase):
    def test_terminal_counts_only_current_completed_tasks_and_keeps_failure_visible(self):
        stream = Terminal()
        tasks = [{'tile_prefix': 'a'}, {'tile_prefix': 'b'}]
        with tempfile.TemporaryDirectory() as root:
            # Old outputs/markers must never count towards a new invocation.
            (Path(root) / '0.done').touch()
            with self.assertRaisesRegex(RuntimeError, 'worker failed'):
                with StageProgress(tasks, 'MPI tiles', root, stream) as progress:
                    progress._count()
                    self.assertEqual(progress.completed, 0)
                    complete_task(progress.tasks[0])
                    progress._count()
                    progress._render(force=True)
                    self.assertEqual(progress.completed, 1)
                    raise RuntimeError('worker failed')
            self.assertFalse(progress.directory.exists())
        self.assertNotIn('_mpi_progress_file', tasks[0])
        self.assertIn('\rMPI tiles: [', stream.getvalue())
        self.assertIn('1/2 | elapsed', stream.getvalue().splitlines()[-1])
        self.assertTrue(stream.getvalue().endswith(' | failed\n'))
        self.assertNotIn(' | complete', stream.getvalue())

    def test_batch_logs_have_throttled_heartbeats_without_terminal_controls(self):
        stream = io.StringIO()
        now = [0.]
        with tempfile.TemporaryDirectory() as root, \
             patch('mesher.mesher_utls.mpi_progress.time.monotonic', side_effect=lambda: now[0]):
            with StageProgress([{}, {}], 'MPI merge stage0 cols', root, stream) as progress:
                complete_task(progress.tasks[0])
                progress._count()
                now[0] = 10
                progress._render()
                self.assertEqual(len(stream.getvalue().splitlines()), 1)
                now[0] = 30
                progress._render()
                self.assertIn('1/2 | elapsed 00:00:30', stream.getvalue())
                # A long-running task still gets a heartbeat with unchanged count.
                now[0] = 60
                progress._render()
                complete_task(progress.tasks[1])
        self.assertEqual(len(stream.getvalue().splitlines()), 4)
        self.assertNotIn('\r', stream.getvalue())
        self.assertNotIn('\x1b', stream.getvalue())
        self.assertTrue(stream.getvalue().endswith('2/2 | elapsed 00:01:00 | complete\n'))

    def test_debug_disables_live_bar_and_separate_stages_have_fresh_markers(self):
        stream = Terminal()
        directories = []
        with tempfile.TemporaryDirectory() as root:
            for _ in range(2):
                with StageProgress([{}], 'MPI tiles', root, stream, debug=True) as progress:
                    directories.append(progress.directory)
                    progress._count()
                    self.assertEqual(progress.completed, 0)
                    complete_task(progress.tasks[0])
        self.assertNotEqual(*directories)
        self.assertNotIn('\r', stream.getvalue())
        self.assertEqual(stream.getvalue().count('1/1'), 2)

    def test_workers_mark_only_successful_tasks(self):
        from mesher.mesher_utls import MPI_mesh_tiles, MPI_merge_tiles
        for module, operation in ((MPI_mesh_tiles, 'mesh_tile'), (MPI_merge_tiles, 'merge_tiles')):
            with self.subTest(worker=module.__name__), tempfile.TemporaryDirectory() as root:
                tasks = [{'_mpi_progress_file': str(Path(root) / f'{i}.done'), 'id': i} for i in range(2)]
                arguments = Path(root) / 'tasks.pickle'
                arguments.write_bytes(cloudpickle.dumps(tasks))

                def work(task):
                    self.assertFalse(Path(task['_mpi_progress_file']).exists())
                    if task['id'] == 1:
                        raise RuntimeError('failed validation')

                with patch.object(module, operation, side_effect=work), \
                     patch.object(module, 'MPI', SimpleNamespace(COMM_WORLD=SimpleNamespace(size=1, rank=0))):
                    with self.assertRaisesRegex(RuntimeError, 'failed validation'):
                        module.main(str(arguments), False)
                self.assertTrue(Path(tasks[0]['_mpi_progress_file']).exists())
                self.assertFalse(Path(tasks[1]['_mpi_progress_file']).exists())
        complete_task({})  # Standalone workers do not need a progress observer.
        with patch('mesher.mesher_utls.mpi_progress.Path.touch', side_effect=OSError('no marker space')):
            complete_task({'_mpi_progress_file': '/unused/0.done'})

    def test_backend_heartbeat_is_debug_only_but_log_is_always_written(self):
        from contextlib import redirect_stdout
        from mesher.mesher_utls.mpi_remesh import _run_backend
        for debug in (False, True):
            with self.subTest(debug=debug), tempfile.TemporaryDirectory() as root:
                stream = io.StringIO()

                def launch(*args, **kwargs):
                    kwargs['stdout'].write('Patch progress: vertices=123\n')
                    kwargs['stdout'].flush()
                    process = MagicMock()
                    process.__enter__.return_value = process
                    process.wait.side_effect = [subprocess.TimeoutExpired('backend', 10), 0]
                    return process

                with redirect_stdout(stream), patch('mesher.mesher_utls.mpi_remesh.subprocess.Popen', side_effect=launch):
                    status, _ = _run_backend(['backend'], Path(root), {'mesher_debug': debug})
                self.assertEqual(status, 0)
                self.assertIn('vertices=123', (Path(root) / 'backend.log').read_text())
                if debug:
                    self.assertIn('Backend active', stream.getvalue())
                    self.assertIn('vertices=123', stream.getvalue())
                else:
                    self.assertEqual(stream.getvalue(), '')

    @unittest.skipUnless(shutil.which('mpirun'), 'OpenMPI is required')
    def test_real_launchers_report_completed_tasks_and_errors(self):
        script = (
            'import json, sys\n'
            'from pathlib import Path\n'
            'from mesher.mesher_cli import _run_mpi_mesh_tasks\n'
            'data = json.loads(Path(sys.argv[1]).read_text())\n'
            '_run_mpi_mesh_tasks(data["tasks"], data["launcher"], 2, "MPI_mesh_tiles.py", "MPI tiles")\n'
        )
        for mode, debug, failure in (('external', False, False), ('spawn', False, False),
                                     ('external', True, False), ('external', False, True)):
            with self.subTest(mode=mode, debug=debug, failure=failure), tempfile.TemporaryDirectory() as root:
                work = Path(root)
                domain = work / 'domain.geojson'
                domain.write_text(json.dumps({'type': 'Polygon', 'coordinates': [
                    [[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]]]}))
                backend = work / 'backend'
                backend.write_text('#!/bin/sh\necho "deliberate backend failure" >&2\nexit 2\n')
                backend.chmod(0o755)
                tasks = [{'tile_prefix': f'tile{i}', 'base_dir': root+'/',
                          'tile_bbox': [0, 0, 5, 5] if failure else [20, 20, 30, 30],
                          'band_width': 1, 'outer_polygon_shp': str(domain), 'constraints': [],
                          'mesher_debug': debug, 'mesher_path': str(backend),
                          'dem_path': str(work/'unused.tif'), 'max_area': 8,
                          'min_area': .25, 'errormetric': 'rmse', 'max_tolerance': 1}
                         for i in range(2)]
                launcher = f'mpirun --oversubscribe -n 2 {shlex.quote(sys.executable)}' if mode == 'external' else None
                args = work / 'input.json'
                args.write_text(json.dumps({'tasks': tasks, 'launcher': launcher}))
                run = subprocess.run([sys.executable, '-c', script, str(args)], cwd=root,
                                     text=True, capture_output=True, timeout=45)
                output = run.stdout + run.stderr
                self.assertIn('MPI tiles: 0/2', output)
                self.assertNotIn('\r', output)
                self.assertFalse(list(work.glob('.mesher-progress-*')))
                if failure:
                    self.assertNotEqual(run.returncode, 0, output)
                    self.assertIn(' | failed', output)
                    self.assertNotIn(' | complete', output)
                    self.assertIn('deliberate backend failure', output)
                    self.assertIn('Reproduction files and full backend log:', output)
                    folders = list(work.glob('cavity-*'))
                    self.assertGreater(len(folders), 0)
                    self.assertIn('deliberate backend failure', (folders[0] / 'backend.log').read_text())
                else:
                    self.assertEqual(run.returncode, 0, output)
                    self.assertIn('MPI tiles: 2/2', output)
                    self.assertIn(' | complete', output)
                    for i in range(2):
                        self.assertTrue((work / f'tile{i}.npz').exists())
                        self.assertTrue(json.loads((work / f'tile{i}.json').read_text())['empty_tile'])
                    self.assertEqual('skipping halo meshing' in output, debug)
