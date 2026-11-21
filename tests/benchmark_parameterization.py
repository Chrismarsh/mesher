"""Compare worker execution with the pre-optimization implementation.

Run from a checkout with the baseline commit available, for example:
PYTHONPATH=src python tests/benchmark_parameterization.py --repeats 3

This is a local single-worker microbenchmark. MPI startup, input preparation,
mesh generation and final exports are excluded from the timed region.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

import cloudpickle
import numpy as np
from osgeo import gdal, osr

from mesher.mesher_utls import MPI_do_parameterize as current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-ref', default='a37a011e378408da43860a691edec01822ef0b80')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--grid-size', type=int, default=100)
    parser.add_argument('--case', choices=('shared-mean', 'batched-stats'), default='shared-mean')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.repeats < 1 or args.grid_size < 1:
        parser.error('repeats and grid-size must be positive')
    if current.MPI.COMM_WORLD.size != 1:
        parser.error('run this single-worker benchmark without mpirun')
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve() if args.output else None
    original_cwd = os.getcwd()
    with tempfile.TemporaryDirectory(prefix='mesher-parameter-benchmark-') as directory:
        work = Path(directory)
        baseline_file = work/'baseline.py'
        baseline_file.write_bytes(subprocess.check_output([
            'git', 'show', f'{args.baseline_ref}:src/mesher/mesher_utls/MPI_do_parameterize.py'], cwd=root))
        spec = importlib.util.spec_from_file_location('baseline_parameterize', baseline_file)
        baseline = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(baseline)
        config = work/'config.py'; config.write_text('mesher_debug=False\n')
        n = args.grid_size
        x, y = np.meshgrid(np.arange(n+1)*2., np.arange(n+1)*2.)
        vertices = np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size)))
        lower = (np.arange(n)[:, None]*(n+1)+np.arange(n)[None, :]).ravel()
        triangles = np.vstack((np.column_stack((lower, lower+1, lower+n+1)),
                               np.column_stack((lower+1, lower+n+2, lower+n+1))))
        mesh = {'mesh': {'vertex': vertices, 'elem': triangles,
                         'nvertex': len(vertices), 'nelem': len(triangles)}}
        srs = osr.SpatialReference(); srs.ImportFromEPSG(32608)
        gt = (0, 1, 0, 2*n, 0, -1)
        raster = work/'values.tif'
        ds = gdal.GetDriverByName('GTiff').Create(str(raster), 2*n, 2*n, 1, gdal.GDT_Float32,
                                                options=['TILED=YES', 'COMPRESS=DEFLATE'])
        ds.SetGeoTransform(gt); ds.SetProjection(srs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(np.random.default_rng(42).random((2*n, 2*n)).astype('f4'))
        ds = None
        methods = ['mean']*3 if args.case == 'shared-mean' else ['mean', 'min', 'max']
        parameters = {f'p{i}': {'filename': [str(raster)], 'method': [method],
                               'classifier': cloudpickle.dumps(lambda value, i=i: value*(i+1))}
                      for i, method in enumerate(methods)}
        parameters['area'] = {}
        old_inputs = [gt, False, mesh, parameters, {}, 2*n, 2*n, srs.ExportToProj4(), True]
        np.save(work/'vertices.npy', vertices)
        np.save(work/'elements.npy', triangles)
        new_inputs = {'gt': gt, 'is_geographic': False, 'parameter_files': parameters,
                      'initial_conditions': {}, 'RasterXSize': 2*n, 'RasterYSize': 2*n,
                      'srs_proj4': srs.ExportToProj4(), 'use_exactextract': True,
                      'nranks': 1, 'output_dir': str(work), 'verts_path': str(work/'vertices.npy'),
                      'elems_path': str(work/'elements.npy')}
        measurements = {'baseline': [], 'optimized': []}
        reference = None
        try:
            os.chdir(work)
            for repeat in range(args.repeats):
                # Alternate execution order to reduce systematic warm-cache bias.
                modes = ['baseline', 'optimized'] if repeat % 2 == 0 else ['optimized', 'baseline']
                for mode in modes:
                    inputs = work/'inputs.pickle'
                    inputs.write_bytes(cloudpickle.dumps(old_inputs if mode == 'baseline' else new_inputs))
                    module = baseline if mode == 'baseline' else current
                    started = time.perf_counter()
                    module.main(str(inputs), False, str(config))
                    measurements[mode].append(time.perf_counter()-started)
                    if mode == 'baseline':
                        result = cloudpickle.loads((work/'pickled_param_args_rets_0.pickle').read_bytes())
                    else:
                        result = cloudpickle.loads((work/'result_0.pickle').read_bytes())['triangles']
                    if reference is None:
                        reference = result
                    if result != reference:
                        raise AssertionError(f'{mode} outputs differ from the reference')
        finally:
            os.chdir(original_cwd)
        medians = {mode: statistics.median(times) for mode, times in measurements.items()}
        report = {'baseline_ref': args.baseline_ref, 'case': args.case, 'triangles': len(triangles),
                  'parameters': 3, 'worker_count': 1, 'repeats': args.repeats,
                  'gdal_version': gdal.VersionInfo(), 'outputs_identical': True,
                  'seconds': measurements, 'median_seconds': medians,
                  'worker_speedup': medians['baseline']/medians['optimized'],
                  'scope': 'Worker execution only; excludes MPI startup, input preparation, meshing and final exports.'}
        text = json.dumps(report, indent=2)+'\n'
        if output:
            output.write_text(text)
        print(text, end='')


if __name__ == '__main__':
    main()
