"""Compact parent-side MPI progress without MPI calls from a background thread."""

import os
from pathlib import Path
import sys
import tempfile
import threading
import time


# Flush per-worker messages only in debug mode so ordinary runs can show a single stage progress
# display.
def debug_log(options, message):
    if options.get('mesher_debug', False):
        print(message, flush=True)


def complete_task(options):
    """Called by a worker only after its mesh and metadata have been written."""
    path = options.get('_mpi_progress_file')
    if path:
        try:
            Path(path).touch()
        except OSError:
            # Losing a progress marker must not turn a valid mesh into a failed
            # job. A successful launcher return still confirms final completion.
            pass


def report_phase(options, message):
    """Let the parent display a long merge's current operation between completions."""
    path = options.get('_mpi_phase_file')
    if path:
        # A single worker owns this file. Replace it atomically so a heartbeat
        # sees a complete message even if it coincides with a phase change.
        target = Path(path)
        scratch = target.with_suffix('.tmp')
        try:
            scratch.write_text(message, encoding='utf-8')
            os.replace(scratch, target)
        except OSError:
            pass  # A display failure must not invalidate a mesh.


class StageProgress:
    """Monitor fresh completion markers on the already-required shared filesystem.

    A single parent thread counts completion markers in one small directory.
    A merge launched with one task can also report its current phase. Neither
    update requires a barrier, a scheduler query or parsing worker logs.
    Temporary arguments/markers are private to this invocation, so old outputs
    cannot inflate completion counts. Backend reproduction files are separate.
    """

    # Copy tasks before assigning run-specific marker paths so the caller can reuse its original
    # task descriptions.
    def __init__(self, tasks, label, directory=None, stream=None, debug=False, show_phase=False):
        self.tasks = [dict(task) for task in tasks]
        self.label = label
        self.root = directory or os.getcwd()
        self.stream = sys.stderr if stream is None else stream
        self.interactive = self.stream.isatty() and not debug
        self.show_phase = show_phase and len(self.tasks) == 1
        self.phase_path = None
        self.completed = 0
        self._stop = threading.Event()
        self._last_line_length = 0
        self._last_update = None

    # Give this launch its own shared marker directory. A resumed stage must count only tasks
    # scheduled in this invocation.
    def __enter__(self):
        self._temporary = tempfile.TemporaryDirectory(prefix='.mesher-progress-', dir=self.root)
        self.directory = Path(self._temporary.name).resolve()
        for i, task in enumerate(self.tasks):
            task['_mpi_progress_file'] = str(self.directory / f'{i}.done')
        if self.show_phase:
            # Thousands of concurrent tiles or pairs already advance the task
            # count. Only a lone merge needs a second signal on shared storage.
            self.phase_path = self.directory / 'merge.phase'
            self.tasks[0]['_mpi_phase_file'] = str(self.phase_path)
        self._started = time.monotonic()
        self._render(force=True)
        self._thread = threading.Thread(target=self._monitor, name='mesher-progress', daemon=True)
        self._thread.start()
        return self

    # Count only completion markers. The same directory also holds launch arguments and can be
    # read while tasks are still running.
    def _count(self):
        with os.scandir(self.directory) as entries:
            self.completed = sum(entry.name.endswith('.done') for entry in entries)

    # A terminal can overwrite its last line. A batch log needs occasional newline-delimited
    # updates that remain readable after a long run.
    def _render(self, force=False, state=None):
        now = time.monotonic()
        # PBS/redirected output gets start/end lines and a 30-second heartbeat;
        # a terminal refreshes one line. Debug output must not fight a live bar.
        interval = 1 if self.interactive else 30
        if not force and self._last_update is not None and now-self._last_update < interval:
            return
        elapsed = int(max(0, now-self._started))
        hours, remainder = divmod(elapsed, 3600)
        minutes, seconds = divmod(remainder, 60)
        duration = f'{hours:02d}:{minutes:02d}:{seconds:02d}'
        total = len(self.tasks)
        bar = ''
        if self.interactive:
            filled = 20*self.completed//total if total else 20
            bar = '[' + '#' * filled + '-' * (20-filled) + '] '
        line = f'{self.label}: {bar}{self.completed}/{total} | elapsed {duration}'
        if state:
            line += f' | {state}'
        elif self.phase_path and self.completed < total:
            try:
                line += ' | ' + self.phase_path.read_text(encoding='utf-8')
            except OSError:
                pass  # The worker may still be starting or publishing its first phase.
        if self.interactive:
            self.stream.write('\r' + line.ljust(self._last_line_length) + ('\n' if state else ''))
        else:
            self.stream.write(line + '\n')
        self.stream.flush()
        self._last_line_length = len(line)
        self._last_update = now

    # Poll from a background thread in the parent process. Workers publish markers without
    # entering extra MPI calls or barriers.
    def _monitor(self):
        while not self._stop.wait(1):
            try:
                self._count()
                self._render()
            except OSError:
                # Progress is optional; a display/filesystem hiccup must not
                # interrupt meshing or replace the actual MPI error.
                return

    # Join the monitor before deleting its directory. On launcher failure, count actual markers
    # so the final line shows how many tasks finished.
    def __exit__(self, exc_type, exc, traceback):
        self._stop.set()
        self._thread.join()
        try:
            if exc_type is None:
                # Successful launcher/Disconnect return confirms every worker
                # finished, even if NFS directory caches lag the last marker.
                self.completed = len(self.tasks)
            else:
                self._count()
            self._render(force=True, state='complete' if exc_type is None else 'failed')
        except OSError:
            pass
        finally:
            self._temporary.cleanup()
