import hashlib
import json
import os
import tempfile
import sys
try:
    from mesher.mesher_utls.bootstrap_utils import check_worker_count, ensure_mesher_on_path
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from mesher.mesher_utls.bootstrap_utils import check_worker_count, ensure_mesher_on_path
ensure_mesher_on_path()
import cloudpickle
from mpi4py import MPI
import numpy as np
import uuid
from osgeo import gdal
import subprocess
from mesher.mesher_utls.raster_utils import fill_nodata

gdal.UseExceptions()  # Enable exception support

def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True

    return False


def projected_path(directory, filename, dataset, warp_options, fill_distance):
    # Include GDAL-reported dependencies (e.g. tiles behind a VRT) without reading raster data.
    dependencies = {}
    pending = list(dataset.GetFileList() or [])
    while pending:
        path = os.path.realpath(pending.pop())
        if path in dependencies:
            continue
        try:
            stat = os.stat(path)
        except OSError:
            dependencies.clear()
            break
        dependencies[path] = (stat.st_size, stat.st_mtime_ns)
        if path.lower().endswith('.vrt') and path != os.path.realpath(filename):
            pending.extend(gdal.Open(path).GetFileList() or [])

    if dependencies:
        recipe = (filename, dependencies, gdal.VersionInfo(), warp_options, fill_distance)
        digest = hashlib.sha256(json.dumps(recipe, sort_keys=True).encode()).hexdigest()[:24]
    else:
        # A GDAL virtual source may have no local file timestamps. Always rebuild it.
        digest = uuid.uuid4().hex
    stem = os.path.splitext(os.path.basename(filename))[0]
    return os.path.join(directory, f'{stem}_{digest}_projected.tif')


def regularize_inputs(args):
    (parameter_dir, data, gdal_prefix, key, pixel_height, pixel_width, srs_out,
     xmax, xmin, ymax, ymin, fill_holes, reuse_projected_parameters) = args


    # there can be multiple files per param output that we use a classifier to merge into one.
    # we need to process each one. Also you can't use a tolerance with the merging classifier as that makes no sense
    if isinstance(data['file'], list) and len(data['file']) == 1:
        print('Error @ ' + key + ': Do not use [ ] for a single file.')
        exit(-1)
    if not 'classifier' in data and isinstance(data['file'], list):
        print(
            'Error @ ' + key + ': If multiple input files are specified for a single parameter a classifer must be provided.')
        exit(-1)
    if 'tolerance' in data and isinstance(data['file'], list):
        print(
            'Error @ ' + key + ': If multiple input files are specified for a single parameter this cannot be used with a tolerance.')
        exit(-1)
    if isinstance(data['file'], list) and not isinstance(data['method'], list):
        print(
            'Error @ ' + key + ': If multiple input files are specified for a single parameter you need to specify each aggregation method.')
        exit(-1)
    if not isinstance(data['file'], list) and isinstance(data['method'], list):
        print(
            'Error @ ' + key + ': If a single input file is specified you cannot specify multiple aggregation methods.')
        exit(-1)
    if not isinstance(data['file'], list):
        data['file'] = [data['file']]
    if not isinstance(data['method'], list):
        data['method'] = [data['method']]
    if len(data['file']) != len(data['method']):
        raise ValueError(f'{key}: Specify one aggregation method per input file.')

    ret_df = dict()
    ret_df['filename'] = []
    ret_df['key'] = key

    fill_distance = int(max(pixel_height, pixel_width) * 5) if fill_holes else None
    for f, resample_method in zip(data['file'], data['method']):
        source = gdal.Open(f)
        if source.GetProjection() == '':
            print("Parameter " + f + " must have spatial reference information.")
            exit(1)

        warp_options = ['-ot', 'Float32', '-overwrite', '-multi', '-dstnodata', '-9999',
                        '-t_srs', srs_out, '-te', *map(str, (xmin, ymin, xmax, ymax)),
                        '-tr', str(pixel_width), str(pixel_height),
                        '-r', 'mode' if resample_method == 'mode' else 'average']
        out_name = projected_path(parameter_dir, f, source, warp_options, fill_distance)
        source = None
        ret_df['filename'].append(out_name)
        if reuse_projected_parameters and os.path.isfile(out_name):
            try:
                if gdal.Open(out_name) is not None:
                    print(f'Reusing projected parameter {key}: {out_name}', flush=True)
                    continue
            except RuntimeError:
                print(f'Rebuilding unreadable projected parameter: {out_name}', flush=True)

        if reuse_projected_parameters == 'resume':
            raise RuntimeError(f'Missing or unreadable projected parameter required for resume: {out_name}')

        # Publish only after warp and optional filling succeed, including when two ranks
        # use the same source and settings for different parameters.
        with tempfile.TemporaryDirectory(dir=parameter_dir) as work_dir:
            tmp = os.path.join(work_dir, 'projected.tif')
            subprocess.check_call([gdal_prefix + 'gdalwarp', *warp_options, f, tmp])
            if fill_holes:
                fill_nodata(tmp, fill_distance)
            os.replace(tmp, out_name)

    return ret_df


def main(pickle_file: str,
         disconnect: bool, expected_workers=None):
    check_worker_count(MPI.COMM_WORLD.size, expected_workers, 'expected GDAL worker count')
    # if called from SLURM, etc, these cli are coming in as strings
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    param_args_split = np.array_split(param_args, MPI.COMM_WORLD.size)

    r = []
    for args in param_args_split[MPI.COMM_WORLD.rank]:
        r.append(regularize_inputs(args))

    # Return this rank's filenames to the parent process.
    with open(os.path.join(os.path.dirname(pickle_file), f'result_{MPI.COMM_WORLD.rank}.pickle'), 'wb') as f:
        cloudpickle.dump(r, f)

    # have been run from the MPI.spawn, so disconnect from parent
    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
gdal.UseExceptions()  # Enable exception support
