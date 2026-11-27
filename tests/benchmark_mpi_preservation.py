"""Compare pairwise remeshing with a historical module and serial patch meshing.

Run with --help. This exercises two tile workers and one pairwise merge, using
identical synthetic raster/quality settings. It is not a cluster scaling test.
Linux /proc supplies sampled process-tree PSS; CPU includes child processes.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
import multiprocessing
import os
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

import numpy as np
from osgeo import gdal
from shapely.geometry import box


# Load a selected remeshing implementation so baseline and candidate run the same input and
# measurement code.
def load_module(filename):
    if not filename:
        from mesher.mesher_utls import mpi_remesh
        return mpi_remesh
    spec = importlib.util.spec_from_file_location('baseline_remesh', filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Mesh one input tile in a child process so the benchmark includes the native backend's work.
def tile_worker(module_path, bounds, options, directory):
    module = load_module(module_path)
    v, t = module.remesh_domain(box(*bounds), [], options, directory)
    return v, t


# Include completed children's CPU time. Timing the Python driver alone would miss much of the
# cost of native refinement.
def cpu_seconds():
    own, children = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


# Generate the tiles and merge them under one configuration so the costs of seam retries and the
# amount of preserved work can be compared.
def run_case(args):
    gdal.UseExceptions()
    with tempfile.TemporaryDirectory(prefix='mesh-benchmark-') as directory:
        work = Path(directory)
        dem = work / 'dem.tif'
        size = args.size
        ds = gdal.GetDriverByName('GTiff').Create(str(dem), size+20, size//2+20, 1, gdal.GDT_Float32)
        ds.SetGeoTransform((-10, 1, 0, size/2+10, 0, -1))
        xx, yy = np.meshgrid(np.arange(size+20)-9.5, size/2+9.5-np.arange(size//2+20))
        ds.GetRasterBand(1).WriteArray(10+80*np.exp(-((xx-size/2)**2+(yy-size/4)**2)/(2*(size/15)**2)))
        ds = None
        baseline = args.mode == 'baseline'
        module_path = args.baseline_module if baseline else None
        module = load_module(module_path)
        options = dict(mesher_path=args.baseline_backend if baseline else args.backend,
                       dem_path=str(dem), max_area=8, min_area=.25, max_tolerance=1,
                       errormetric='rmse', remesh_band_width=2, remesh_max_attempts=12,
                       remesh_max_vertices=500000)
        domain = box(0, 0, size, size/2)
        initial_cpu = cpu_seconds()
        started = time.perf_counter()
        if args.mode == 'serial':
            result = module.remesh_domain(domain, [], options, work)
            v, t = result
            stats = getattr(result, 'stats', {})
            tile_seconds = 0
        else:
            with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context('spawn')) as pool:
                futures = [pool.submit(tile_worker, module_path, bounds, options, directory) for bounds in
                           ((0, 0, size/2+3, size/2), (size/2-3, 0, size, size/2))]
                a, b = [future.result() for future in futures]
            tile_seconds = time.perf_counter()-started
            calls = []
            # Record the area requested and the number of triangles generated for each real
            # patch solve so repeat work is visible in the benchmark.
            def measured_patch(cavity, constraints, options, directory, fixed):
                call = {'area': cavity.area}
                calls.append(call)
                result = module.remesh_domain(cavity, constraints, options, directory, fixed)
                call['triangles'] = len(result[1])
                return result
            v, t, stats = module.merge_cavity(*a, [0, 0, size/2, size/2], *b, [size/2, 0, size, size/2],
                                             domain, [], options, directory, measured_patch)
            stats['cavity_calls'] = calls
            stats['total_attempted_cavity_area'] = sum(c['area'] for c in calls)
            stats['generated_patch_triangles'] = sum(c.get('triangles', 0) for c in calls)
        elapsed = time.perf_counter()-started
        used_cpu = cpu_seconds()-initial_cpu
        # Independent acceptance checks, excluded from measured generation time.
        from mesher.mesher_utls.mpi_remesh import validate_patch, _stitch
        validate_patch(v, t, domain)
        _stitch(v, t, np.empty((0, 3)), np.empty((0, 3), dtype=np.int64), domain)
        Path(args.output).write_text(json.dumps(dict(mode=args.mode, size=size, wall_seconds=elapsed,
            tile_seconds=tile_seconds, merge_seconds=elapsed-tile_seconds if args.mode != 'serial' else 0,
            cpu_seconds=used_cpu,
            triangles=len(t), vertices=len(v), stats=stats), indent=2))


# Sum proportional set sizes across the process tree so shared mappings are apportioned across
# workers in the memory estimate.
def tree_pss(root):
    parents = {}
    for directory in Path('/proc').glob('[0-9]*'):
        try:
            fields = (directory / 'stat').read_text().rsplit(')', 1)[1].split()
            parents[int(directory.name)] = int(fields[1])
        except (OSError, ValueError, IndexError):
            pass
    selected = {root}
    while True:
        children = {pid for pid, parent in parents.items() if parent in selected}
        if children.issubset(selected): break
        selected |= children
    total = 0
    for pid in selected:
        try:
            for line in Path(f'/proc/{pid}/smaps_rollup').read_text().splitlines():
                if line.startswith('Pss:'): total += int(line.split()[1])*1024
        except (OSError, ValueError):
            pass
    return total


# Sample memory from outside the measured run. The driver remains able to observe workers while
# their native backends are busy.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', required=True)
    parser.add_argument('--baseline-backend', required=True)
    parser.add_argument('--baseline-module', required=True, help='Historical mpi_remesh.py exported with git show')
    parser.add_argument('--size', type=int, default=400)
    parser.add_argument('--output', required=True, help='JSON result path; logs saved alongside it')
    parser.add_argument('--mode', choices=['baseline', 'preserved', 'serial'], help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.mode:
        run_case(args)
        return
    results = []
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    for mode in ('baseline', 'preserved', 'serial'):
        result_file = output.with_name(output.stem+'-'+mode+'.json')
        command = [sys.executable, str(Path(__file__).resolve()), '--backend', args.backend,
                   '--baseline-backend', args.baseline_backend, '--baseline-module', args.baseline_module,
                   '--size', str(args.size), '--output', str(result_file), '--mode', mode]
        with output.with_name(output.stem+'-'+mode+'.log').open('w') as log:
            with subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT) as process:
                peak = 0
                while process.poll() is None:
                    peak = max(peak, tree_pss(process.pid))
                    time.sleep(.2)
                if process.returncode:
                    raise RuntimeError(f'{mode} failed; see its log')
        result = json.loads(result_file.read_text())
        result['sampled_peak_process_tree_pss_bytes'] = peak
        results.append(result)
        print(f'{mode}: {result["wall_seconds"]:.3f}s, {result["triangles"]} triangles', flush=True)
    output.write_text(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
