"""Measure export heap growth (file-cache RSS is intentionally reported separately).

Run with PYTHONPATH=src python tests/benchmark_streaming.py --faces 100000 1000000.
Uses repeated valid cells to isolate storage/export cost, not topology/meshing.
"""
import argparse
import json
import resource
import tempfile
import time
import tracemalloc
from pathlib import Path

import numpy as np
from mesher.mesh_store import batches, disk_array, dump_json, write_vtu


def measure(count):
    with tempfile.TemporaryDirectory(prefix='mesher-stream-benchmark-') as work:
        root = Path(work)
        vertices = disk_array(root/'vertices.npy', (3, 3))
        vertices[:] = [[0, 0, 1], [1, 0, 2], [0, 1, 3]]
        triangles = disk_array(root/'triangles.npy', (count, 3), np.int64)
        values = disk_array(root/'values.npy', (count,))
        for lo, hi in batches(count):
            triangles[lo:hi] = [0, 1, 2]
            values[lo:hi] = np.arange(lo, hi)
        mesh = {'mesh': dict(vertex=vertices, elem=triangles, nelem=count, nvertex=3, proj4='')}
        tracemalloc.start()
        started = time.perf_counter()
        with open(root/'mesh.json', 'w') as stream:
            dump_json(mesh, stream)
        with open(root/'params.json', 'w') as stream:
            dump_json({'value': values}, stream)
        write_vtu(root/'mesh.vtu', mesh, {'value': values}, {})
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return dict(faces=count, seconds=time.perf_counter()-started,
                    python_numpy_peak_bytes=peak,
                    process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--faces', type=int, nargs='+', default=[100000, 1000000])
    args = parser.parse_args()
    for faces in args.faces:
        print(json.dumps(measure(faces)), flush=True)
