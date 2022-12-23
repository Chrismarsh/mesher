import os
import sys
import cloudpickle
from mpi4py import MPI
import numpy as np
import uuid
from osgeo import gdal
import subprocess

def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True

    return False

def regularize_inputs(args):
    base_dir, data, exec_str, gdal_prefix, key, pixel_height, pixel_width, srs_out, xmax, xmin, ymax, ymin, fill_holes = args

    # make a copy as this exec string is reused for inital conditions
    # and we don't want the changes made here to impact it
    estr = exec_str
    if data['method'] == 'mode':
        estr = exec_str + 'mode'
    else:
        estr = exec_str + 'average'

    do_cell_resize = True
    estr = estr + ' -tr %s %s'
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

    ret_df = dict()
    ret_df['filename'] = []
    ret_df['key'] = key

    for f in data['file']:
        # we need to handle a path being passed in
        output_param_fname = os.path.basename(f)
        output_param_fname = os.path.splitext(output_param_fname)[0]

        if gdal.Open(f).GetProjection() == '':
            print("Parameter " + f + " must have spatial reference information.")
            exit(1)

        mangle = uuid.uuid4().hex[:8]

        out_name = base_dir + output_param_fname + '_' + mangle + '_projected.tif'

        # force all the parameter files to have the same extent as the input DEM
        subprocess.check_call([estr % (gdal_prefix,
                                       f, out_name, srs_out, xmin, ymin, xmax, ymax, pixel_width,
                                       pixel_height)], shell=True)
        if fill_holes:
            dist = max([pixel_height,pixel_width]) * 5
            exec_str = '%sgdal_fillnodata.py %s -md %d' % (gdal_prefix, out_name, dist)
            subprocess.check_call([exec_str], shell=True)

        ret_df['filename'].append(out_name)

    return ret_df


def main(pickle_file: str,
         disconnect: bool):
    # if called from SLURM, etc, these cli are coming in as strings
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    param_args_split = np.array_split(param_args, MPI.COMM_WORLD.size)

    r = []
    for args in param_args_split[MPI.COMM_WORLD.rank]:
        r.append(regularize_inputs(args))

    # there is no way to return the uuid mangled filename + param name  so save it to a pickly
    with open(f'pickled_param_args_rets_{MPI.COMM_WORLD.rank}.pickle', 'wb') as f:
        cloudpickle.dump(r,f)

    # have been run from the MPI.spawn, so disconnect from parent
    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
