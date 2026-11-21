"""Compact parent-side MPI progress without MPI calls from a background thread."""

import os
from pathlib import Path
import sys
import tempfile
import threading
import time


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


class StageProgress:
    """Monitor fresh completion markers on the already-required shared filesystem.

    A single parent thread reads one small directory per second. Workers write
    once per completed task, with no barriers, scheduler queries or log parsing.
    Temporary arguments/markers are private to this invocation, so old outputs
    cannot inflate completion counts. Backend reproduction files are separate.
    """

    def __init__(self, tasks, label, directory=None, stream=None, debug=False):
        self.tasks = [dict(task) for task in tasks]
        self.label = label
        self.root = directory or os.getcwd()
        self.stream = sys.stderr if stream is None else stream
        self.interactive = self.stream.isatty() and not debug
        self.completed = 0
        self._stop = threading.Event()
        self._last_line_length = 0
        self._last_update = None

    def __enter__(self):
        self._temporary = tempfile.TemporaryDirectory(prefix='.mesher-progress-', dir=self.root)
        self.directory = Path(self._temporary.name).resolve()
        for i, task in enumerate(self.tasks):
            task['_mpi_progress_file'] = str(self.directory / f'{i}.done')
        self._started = time.monotonic()
        self._render(force=True)
        self._thread = threading.Thread(target=self._monitor, name='mesher-progress', daemon=True)
        self._thread.start()
        return self

    def _count(self):
        with os.scandir(self.directory) as entries:
            self.completed = sum(entry.name.endswith('.done') for entry in entries)

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
        if self.interactive:
            self.stream.write('\r' + line.ljust(self._last_line_length) + ('\n' if state else ''))
        else:
            self.stream.write(line + '\n')
        self.stream.flush()
        self._last_line_length = len(line)
        self._last_update = now

    def _monitor(self):
        while not self._stop.wait(1):
            try:
                self._count()
                self._render()
            except OSError:
                # Progress is optional; a display/filesystem hiccup must not
                # interrupt meshing or replace the actual MPI error.
                return

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
