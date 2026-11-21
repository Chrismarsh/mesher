#!/usr/bin/env python
# Mesher
# Copyright (C) 2017 Christopher Marsh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from osgeo import gdal, ogr, osr
import subprocess
import re
import json
import os
import shlex
import tempfile
import numpy as np
import math

import sys
import shutil
import importlib
import warnings
import cloudpickle

import time
from mpi4py import MPI
import glob



from . import mesher_utls
from .mesh_store import (BATCH_SIZE, batches, disk_array, copy_array, open_mesh_arrays,
                         check_boundary_geometry, dump_json, write_vtu)
from .mesh_topology import (MeshTopologyError, integer_rows, build_neighbours,
                            boundary_subdivisions)
from .mesher_utls.raster_utils import fill_nodata
from .mesher_utls.mpi_progress import StageProgress

gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()  # Enable exception support
osr.UseExceptions()  # Enable exception support

# Version of the mesh that is produced
MESH_MAJOR = "1"
MESH_MINOR = "0"
MESH_PATCH = "0"


def main():
    # load user configurable parameters here
    # Check user defined configuration file

    if len(sys.argv) == 1:
        print('ERROR: meshgen.py requires one argument [configuration file] (i.e. meshgen.py Bow.py)')
        return

    run_started = time.perf_counter()
    # Get name of configuration file/module
    configfile = sys.argv[-1]

    # Read configuration
    X, bufferDist, clip_to_shp, constraints, dem_filename, do_smoothing, errormetric, extent, fill_holes, \
        initial_conditions, lloyd_itr, max_area, max_smooth_iter, max_tolerance, mesher_path, no_simplify_buffer,\
        nworkers, nworkers_gdal, output_write_shp, output_write_vtu, parameter_files, reuse_mesh, scaling_factor, simplify,\
        simplify_tol, use_input_prj, user_no_weights, user_output_dir, verbose, weight_threshold, wkt_out, \
        MPI_exec_str, MPI_nworkers, mpi_mesh, \
        mesher_debug, use_exactextract, write_param_per_file = read_config(configfile)


    ########################################################

    mpi_remesh_options = mpi_remesh_options_from_config(X) if mpi_mesh else None
    num_tiles = getattr(X, 'mpi_num_tiles', None) if mpi_mesh else None
    if num_tiles is not None:
        _positive_integer(num_tiles, 'mpi_num_tiles')
    parameter_batch_size = _positive_integer(getattr(X, 'mpi_parameter_batch_size', BATCH_SIZE),
                                             'mpi_parameter_batch_size')
    gdal_prefix = find_gdal_prefix()

    base_name = os.path.basename(dem_filename)
    base_name = os.path.splitext(base_name)[0]

    base_dir = user_output_dir + base_name + os.path.sep
    parameter_dir = os.path.join(base_dir, 'projected_parameters')
    reuse_projected_parameters = getattr(X, 'reuse_projected_parameters', False)

    if not os.path.isdir(base_dir) and reuse_mesh:
        raise Exception(f"Using reuse_mesh=True, however the directory {base_dir} does not exist")

    # A fresh mesh needs fresh intermediates, even when the parameter rasters are reused.
    if os.path.isdir(base_dir) and not reuse_mesh:
        if os.path.islink(os.path.normpath(base_dir)):
            raise OSError(f'Cannot clean a symlinked mesh directory: {base_dir}')
        for entry in os.scandir(base_dir):
            if entry.is_dir(follow_symlinks=False):
                if reuse_projected_parameters and entry.path == parameter_dir:
                    continue
                shutil.rmtree(entry.path)
            else:
                os.remove(entry.path)

    os.makedirs(base_dir, exist_ok=True)

    # we want to reuse an already generated mesh, but we will need to clean up the shp file as gdal won't overwrite
    # an existing one
    if reuse_mesh:
        try:
            os.remove(base_dir + base_name + '_USM.shp')
        except OSError:
            pass  # If the folder doesn't exist, don't worry

    try:
        src_ds = gdal.Open(dem_filename)
    except RuntimeError as e:
        print('Unable to open file ' + dem_filename)
        raise e

    if src_ds.GetProjection() == '':
        raise RuntimeError("Input DEM must have spatial reference information.")

    # if we are using the input DEM's projection, grab it here from the source DEM
    if use_input_prj:
        wkt_out = src_ds.GetProjection()

    srs_out = osr.SpatialReference()

    srs_out.ImportFromWkt(wkt_out)

    if do_smoothing:
        output_file_name = '_projected_0.tif'
    else:
        output_file_name = '_projected.tif'

    ext_str = ''

    # clip to a user-specified extent
    if extent is not None:
        wkt = src_ds.GetProjection()
        srs = osr.SpatialReference()

        srs.ImportFromWkt(wkt)

        ext_str = ' -te %s %s %s %s -te_srs \"%s\" ' % (extent[0], extent[1], extent[2], extent[3], srs.ExportToProj4())

    cutline_cmd = ''
    if clip_to_shp:
        cutline_cmd = f' -cutline {clip_to_shp} -crop_to_cutline '

    e = 'GDAL_CACHEMAX=\"5%%\" %sgdalwarp %s %s -ot Float32 -multi -overwrite -dstnodata -9999 %s -t_srs \"%s\"' + ext_str
    subprocess.check_call([e % (gdal_prefix,
                                dem_filename, base_dir + base_name + output_file_name, cutline_cmd, srs_out.ExportToProj4())],
                          shell=True)

    # fill NoData holes in the source DEM
    # default is 5 * resolution
    if fill_holes:
        gt = src_ds.GetGeoTransform()

        pixel_width = gt[1]
        pixel_height = gt[5]

        dist = max([pixel_height, pixel_width]) * 5
        fill_nodata(base_dir + base_name + output_file_name, int(dist))

    src_ds = gdal.Open(base_dir + base_name + output_file_name)

    if src_ds is None:
        print('Unable to open %s' % dem_filename)
        exit(1)

    #obtain new values for this after all the projection, etc
    gt = src_ds.GetGeoTransform()

    pixel_width = gt[1]
    pixel_height = gt[5]

    if hasattr(X, 'min_area'):
        min_area = X.min_area
    else:
        # if the user doesn't specify, then limit to the underlying resolution. No point going past this!
        min_area = abs(
            pixel_width * pixel_height)


    # DEM smoothing passes
    if do_smoothing:
        for itr in range(max_smooth_iter):

            in_name = base_dir + base_name + '_projected_%d.tif' % itr
            out_name = base_dir + base_name + '_projected_%d.tif' % (itr + 1)

            if itr + 1 == max_smooth_iter:  # last iteration, change output
                out_name = base_dir + base_name + '_projected.tif'

            subprocess.check_call(
                [
                    'GDAL_CACHEMAX=\"5%%\" %sgdalwarp %s %s -ot Float32 -multi  -overwrite -dstnodata -9999 -r cubicspline -tr %s %s' % (
                        gdal_prefix,
                        in_name, out_name, abs(pixel_width) / scaling_factor,
                        abs(pixel_height) / scaling_factor)], shell=True)

            scaling_factor *= (itr + 1)

    # now, reopen the file
    src_ds = gdal.Open(base_dir + base_name + '_projected.tif')

    # create the spatial reference from the raster dataset
    wkt = src_ds.GetProjection()
    srs = osr.SpatialReference()

    srs.ImportFromWkt(wkt)

    is_geographic = srs.IsGeographic()

    if is_geographic and simplify:
        print('Error: Simplify=True and a geographic output CRS is not currently supported.')
        exit(1)

    gt = src_ds.GetGeoTransform()
    # x,y origin
    xmin = gt[0]
    ymax = gt[3]

    pixel_width = gt[1]
    pixel_height = gt[5]

    xmax = xmin + pixel_width * src_ds.RasterXSize
    ymin = ymax + pixel_height * src_ds.RasterYSize  # pixel_height is negative

    total_weights_param, use_weights_param = regularize_inputs(parameter_dir, gdal_prefix, parameter_files,
                                                               pixel_height, pixel_width,
                                                               srs_out, xmax, xmin, ymax, ymin, fill_holes,
                                                               MPI_exec_str, MPI_nworkers, reuse_projected_parameters)

    total_weights_ic, use_weights_ic = 0, False

    use_weights = use_weights_param or use_weights_ic

    topo_weight = 1 - (total_weights_param + total_weights_ic)

    # over ride with user settings
    if user_no_weights:
        use_weights = False

    if use_weights and topo_weight < 0:
        raise RuntimeError("Parameter weights must equal 1")

    # the empty area to trigger area computations later
    # needs to happen here after the regularization as this doesn't have anything to regularize
    parameter_files['area'] = {}

    plgs_shp = base_name + '.shp'

    dem = src_ds.GetRasterBand(1)

    # Create a mask raster that has a uniform value for all cells read all elevation data. Might be worth changing to
    # the for-loop approach we use below so we don't have to read in all into ram. some raster
    Z = dem.ReadAsArray()

    # set all the non-nodata values to our mask value
    Z[Z != dem.GetNoDataValue()] = 1

    # save to file
    (x, y) = Z.shape
    driver = gdal.GetDriverByName("GTiff")
    dst_datatype = gdal.GDT_Int16  # gdal.GDT_Float32   #<--- this is fine as we don't actually store elevation data.
    tmp_raster = base_dir + 'mask_' + base_name + '.tif'
    dst_ds = driver.Create(tmp_raster, y, x, 1, dst_datatype)
    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjection())
    dst_ds.GetRasterBand(1).SetNoDataValue(dem.GetNoDataValue())
    dst_ds.GetRasterBand(1).WriteArray(Z)
    dst_ds.FlushCache()  # super key to get this to disk
    # noinspection PyUnusedLocal
    dst_ds = None  # close file

    # raster -> polygon
    # subprocess.check_call(
    #     ['%sgdal_polygonize.py %s -b 1 -mask %s -f "ESRI Shapefile" %s' % (gdal_prefix, tmp_raster, tmp_raster,
    #                                                                        base_dir +
    #                                                                        plgs_shp)], shell=True)

    gdal_polygonize(tmp_raster, tmp_raster, base_dir + plgs_shp)

    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(base_dir + plgs_shp, 1)

    # If the input has multiple polygon regions, find the largest polygon and use that as the meshing domain
    layer = dataSource.GetLayer()
    max_geom_area = -1
    max_feature_ID = None
    for feature in layer:
        geom = feature.GetGeometryRef()
        area = geom.GetArea()
        # print 'FID = ' + str(feature.GetFID()) + ' area = ' + str(area)
        if area > max_geom_area:
            max_feature_ID = feature.GetFID()
            max_geom_area = area

    print('Using FID = ' + str(max_feature_ID) + " as the largest continuous area.")
    feats = np.arange(0, layer.GetFeatureCount())
    for f in feats:
        if f != max_feature_ID:
            layer.DeleteFeature(f)
    layer.SyncToDisk()
    dataSource.ExecuteSQL("REPACK " + layer.GetName())
    dataSource.Destroy()

    # allow us to be able to use this directly in the exec_string below if we don't simplify
    outputBufferfn = base_dir + plgs_shp

    # simplify the outter domain constraint and contract by X m, which allows more flexibility in fitting larger
    # triangles to the outter domain.
    if simplify and not no_simplify_buffer:
        # buffering code from http://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html?highlight=buffer
        inputfn = base_dir + plgs_shp
        outputBufferfn = base_dir + "buffered_" + plgs_shp

        print('Simplifying extents by buffer distance = ' + str(bufferDist))
        inputds = ogr.Open(inputfn)
        inputlyr = inputds.GetLayer()

        shpdriver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(outputBufferfn):
            shpdriver.DeleteDataSource(outputBufferfn)
        outputBufferds = shpdriver.CreateDataSource(outputBufferfn)
        bufferlyr = outputBufferds.CreateLayer(outputBufferfn, srs_out, geom_type=ogr.wkbPolygon)
        featureDefn = bufferlyr.GetLayerDefn()

        for feature in inputlyr:
            raw_ingeom = feature.GetGeometryRef()
            # collapse the raster-resolution staircase (e.g. a 500km square on a 15m grid is
            # ~130k collinear vertices) so the buffer below is cheap. Accuracy is bounded by
            # simplify_tol, which is the same tolerance the linestring conversion uses later.
            # NOTE: SimplifyPreserveTopology() returns a NEW geometry, it does not modify the
            # input in place, so the returned geometry must be used.
            ingeom = raw_ingeom
            simplified = None
            try:
                simplified = raw_ingeom.SimplifyPreserveTopology(simplify_tol)
            except Exception as e:
                print('Warning: SimplifyPreserveTopology failed (%s). Using original boundary.' % e)
            if simplified is not None and not simplified.IsEmpty():
                ingeom = simplified
            elif simplified is not None:
                print('Warning: simplified boundary is empty. Using original boundary.')
            geomBuffer = ingeom.Buffer(bufferDist)
            print('Outer boundary vertices: %d raw -> %d simplified -> %d buffered' %
                  (_geom_point_count(raw_ingeom), _geom_point_count(ingeom),
                   _geom_point_count(geomBuffer)))

            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(geomBuffer)
            bufferlyr.CreateFeature(outFeature)
            outFeature = None
        # close the files
        inputds = None
        outputBufferds = None

    # this needs to be done after we optionally simplify the outer domain. If we do, we need to ensure any interior
    # constraints are clipped to this new domain as the triangulation does not like intersecting constraints.
    for key, data in constraints.items():
        # we need to handle a path being passed in
        output_constraint_fname = os.path.basename(data['file'])
        output_constraint_fname = os.path.splitext(output_constraint_fname)[0]

        df = ogr.Open(data['file'])
        if df is None:
            raise RuntimeError("Unable to open constraint " + data['file'])

        if df.GetLayer(0).GetSpatialRef() is None:
            raise RuntimeError("Constraint " + data['file'] + " must have spatial reference information.")
        df = None

        outname = base_dir + 'constraint_' + output_constraint_fname

        # force all the constraints to have the same extent as the input DEM
        exec_string = '%sogr2ogr -overwrite %s %s  -t_srs \"%s\"' % (
            gdal_prefix, outname + '.shp', data['file'], wkt_out.replace('"', '\\"'))

        if 'simplify' in data:
            exec_string = exec_string + ' -simplify ' + str(
                data['simplify'])  # because of ogr2ogr, the simplification is done in the units of the original data

        subprocess.check_call(exec_string, shell=True)

        # if we simplified with a buffer, clip the constraint shp to that extent. Do this after we project it above.
        if simplify and not no_simplify_buffer:
            clip_outname = base_dir + 'clip_constraint_' + output_constraint_fname
            exec_string = '%sogr2ogr -f "ESRI Shapefile" -clipsrc %s %s %s' % (
                gdal_prefix, outputBufferfn, clip_outname + '.shp', outname + '.shp')  # clip src, output, input
            subprocess.check_call(exec_string, shell=True)

            # update outname to be the clipped one so we can use it below
            outname = clip_outname

        # convert to geoJSON because it's easy to parse
        # ensure it's all line strings and explode the multilines into linestrings
        subprocess.check_call(['%sogr2ogr -f GeoJSON   -nlt LINESTRING -explodecollections  %s %s' % (
            gdal_prefix, outname + '.geojson', outname + '.shp')], shell=True)

        constraints[key]['filename'] = outname + '.geojson'

        # since all the project has already happened in the ogr2ogr step, we can just read it in now
        with open(outname + '.geojson') as f:
            constraints[key]['file'] = json.load(f)

    print('Converting polygon to linestring')
    exec_string = '%sogr2ogr -overwrite %s %s  -nlt LINESTRING' % (
        gdal_prefix, base_dir + 'line_' + plgs_shp, outputBufferfn)

    if simplify:
        exec_string = exec_string + ' -simplify ' + str(simplify_tol)

    subprocess.check_call(exec_string, shell=True)

    # convert to geoJSON because it's easy to parse
    poly_plgs = base_name + '.geojson'
    subprocess.check_call(['%sogr2ogr -f GeoJSON %s %s' % (gdal_prefix, base_dir +
                                                           poly_plgs, base_dir + 'line_' + plgs_shp)], shell=True)

    with open(base_dir + poly_plgs) as f:
        plgs = json.load(f)

    # look through all the features and find the biggest
    idx = -1
    i = 0
    cmax = -1
    l = 0

    for features in plgs['features']:
        if features['geometry'] is None:
            print(
                'Error: the choice of simplify buffer and tolerance has resulted in oversimplifying the domain. '
                'Please choose a tighter tolerance')
            exit(1)

        if features['geometry']['type'] == 'LineString':
            l = len(features['geometry']['coordinates'])
        elif features['geometry']['type'] == 'MultiLineString':
            coords = []

            idx_ml = 0
            len_ml = -1
            j = 0
            # find the largest of the multi lines, ignore holes!
            for lines in features['geometry']['coordinates']:
                l = len(lines)
                if l > len_ml:
                    len_ml = l
                    idx_ml = j
                j += 1

            for l in features['geometry']['coordinates'][idx_ml]:
                coords.append(l)

            features['geometry']['type'] = 'LineString'
            features['geometry']['coordinates'] = coords
            l = len(coords)

        if l > cmax:
            cmax = l
            idx = i
        i += 1

    # assuming just the biggest feature is what we want. Need to add in more to support rivers and lakes
    if plgs['features'][idx]['geometry']['type'] != 'LineString':
        raise RuntimeError('Not linestring')

    coords = plgs['features'][idx]['geometry']['coordinates']

    # We can't just insert this into the global PLGS datastructure as we need that to define the convex hull so merge
    # all the constraints into 1 geojson file so we can just load in the main cpp mesher code to define the interior
    # PLGS TODO: in the future this and the convex hull PLGS should all be in 1 PLGS and the Triangle compatibility
    #  should be removed
    interior_PLGS = {
        "type": "FeatureCollection",
        "name": "interior_PLGS",
        "features": []}

    for key, data in constraints.items():  # over each constraint
        for feat in data['file']['features']:  # over the features present in each constraint

            # round down to 5 decimal places. This helps ensure points that differ by round-off levels aren't treated
            # as seperate points
            # Since we are in projected units (metres??) this should almost certainly be fine...
            # todo: check this is robust
            feat['geometry']['coordinates'] = [[round(item, 5) for item in sublist] for sublist in feat['geometry']['coordinates']]
            interior_PLGS['features'].append(feat)

    with open(base_dir + 'interior_PLGS.geojson', 'w') as fp:
        json.dump(interior_PLGS, fp)

    # Create the outer PLGS to constrain the triangulation
    poly_file = 'PLGS' + base_name + '.poly'
    with open(base_dir + poly_file, 'w') as f:
        header = '%d 2 0 0\n' % (len(coords))
        f.write(header)
        vert = 1
        for c in coords:
            f.write('%d %17.11f %17.11f\n' % (vert, c[0], c[1]))
            vert = vert + 1

        f.write('\n')
        header = '%d 0\n' % (len(coords))
        f.write(header)

        for i in range(len(coords)):
            if i + 1 == len(coords):  # last time has to loop back to start
                f.write('%d %d %d\n' % (i + 1, i + 1, 1))
            else:
                f.write('%d %d %d\n' % (i + 1, i + 1, i + 2))

        f.write('0\n')

    print(f'Preprocessing took {time.perf_counter() - run_started:.2f} s')
    meshing_started = time.perf_counter()
    final_npz_path = None
    if mpi_mesh and reuse_mesh:
        final_npz_path = compute_mpi_final_npz_path(base_dir, base_name, MPI_nworkers)
        if not os.path.exists(final_npz_path):
            raise RuntimeError(
                f'reuse_mesh=True and mpi_mesh=True, but stitched mesh not found at {final_npz_path}. '
                'Ensure the previous MPI mesh run completed; legacy outputs without a layout manifest require the original MPI_nworkers.'
            )
        validate_mpi_mesh_metadata(final_npz_path)

    # if we aren't reusing the mesh, generate a new one
    if not reuse_mesh:
        print(f'mpi_mesh={mpi_mesh}')
        if mpi_mesh:
            final_npz_path = run_mpi_meshing(base_dir, base_name, xmin, ymin, xmax, ymax,
                                             mesher_path, outputBufferfn, constraints, parameter_files,
                                             initial_conditions, max_area, min_area, max_tolerance,
                                             errormetric, use_weights, topo_weight,
                                             weight_threshold, is_geographic, MPI_exec_str, MPI_nworkers,
                                             mesher_debug=mesher_debug, remesh_options=mpi_remesh_options,
                                             num_tiles=num_tiles)
        else:
            execstr = '%s --poly-file %s --tolerance %s --raster %s --area %s --min-area %s --error-metric %s --lloyd %d --interior-plgs-file %s' % \
                      (mesher_path,
                       base_dir + poly_file,
                       max_tolerance,
                       base_dir + base_name + '_projected.tif',
                       max_area,
                       min_area,
                       errormetric,
                       lloyd_itr,
                       base_dir + 'interior_PLGS.geojson'
                       )

            if is_geographic:
                execstr += ' --is-geographic true'

            if use_weights:
                execstr += ' --weight %s' % topo_weight
                execstr += ' --weight-threshold %s' % weight_threshold
            if mesher_debug:
                execstr += ' --debug true'

            for key, data in parameter_files.items():
                if 'tolerance' in data:
                    if data['method'] == 'mode':
                        execstr += ' --category-raster %s --category-frac %s' % (data['filename'][0], data[
                            'tolerance'])  # we need [0] on raster as it's been made iterable by this point
                    else:
                        execstr += ' --raster %s --tolerance %s' % (data['filename'][0], data['tolerance'])
                if use_weights and 'weight' in data:
                    execstr += ' --weight %s' % data['weight']

            for key, data in initial_conditions.items():
                if 'tolerance' in data:
                    if data['method'] == 'mode':
                        execstr += ' --category-raster %s --category-frac %s' % (data['filename'], data['tolerance'])
                    else:
                        execstr += ' --raster %s --tolerance %s' % (data['filename'][0], data['tolerance'])
                if use_weights and 'weight' in data:
                    execstr += ' --weight %s' % data['weight']

            print(execstr)
            subprocess.check_call(execstr, shell=True)

    mesh_action = 'Mesh reuse' if reuse_mesh else 'Meshing'
    print(f'{mesh_action} took {time.perf_counter() - meshing_started:.2f} s')
    mesh_preparation_started = time.perf_counter()

    # some paramters we want to use to constrain the mesh but don't actually want in the output. This let's use
    # remove them
    keys_to_drop = []
    for key, data in parameter_files.items():
        if 'drop' in data and data['drop'] == True:
            keys_to_drop.append(key)
    for key in keys_to_drop:
        print('Dropping parameter ' + key + ' and not writing to mesh due to user request')
        del parameter_files[key]

    # The preparation store lives through export; workers map the same files.
    with tempfile.TemporaryDirectory(prefix='mesher-prepared-', dir=base_dir) as preparation_dir:
        if mpi_mesh:
            if final_npz_path is None:
                raise RuntimeError('mpi_mesh=True but no stitched mesh was produced')
            with open_mesh_arrays(final_npz_path, preparation_dir) as (verts, tris):
                mesh = build_mesh_from_arrays(verts, tris, src_ds, dem, srs, is_geographic, verbose,
                                              directory=preparation_dir)
        else:
            mesh = load_mesh_from_mesher_files(base_dir, base_name, src_ds, dem, srs,
                                               is_geographic, verbose, directory=preparation_dir)

        print(f'Mesh preparation took {time.perf_counter() - mesh_preparation_started:.2f} s')
        ics = {key: [] for key in initial_conditions}

        gt = src_ds.GetGeoTransform()


        print('Computing parameters and initial conditions')

        # prepare the meta data bundle so the parallel ranks can produce their subset to work on in parallel
        start_time = time.perf_counter()

        parameter_args = {
            'gt': gt,
            'is_geographic': is_geographic,
            'parameter_files': parameter_files,
            'initial_conditions': initial_conditions,
            'RasterXSize': src_ds.RasterXSize,
            'RasterYSize': src_ds.RasterYSize,
            'srs_proj4': srs_out.ExportToProj4(),
            'use_exactextract': use_exactextract,
            'batch_size': parameter_batch_size,
        }
        params = run_parameterization(mesh, parameter_args, base_dir,
                                      MPI_exec_str, MPI_nworkers, configfile,
                                      output_dir=preparation_dir)

        for key, data in parameter_files.items():
            parameter_files[key]['file'] = None

        for key, data in initial_conditions.items():
            initial_conditions[key]['file'] = None

        mesh['mesh']['version'] = f'{MESH_MAJOR}.{MESH_MINOR}.{MESH_PATCH}'

        print(f'Parameterization took {time.perf_counter() - start_time:.2f} s')

        export_started = time.perf_counter()
        output_name = base_name + '_mpi' if mpi_mesh else base_name

        if output_write_vtu:
            print('Writing vtu file...')
            write_vtu(base_dir + output_name + '.vtu', mesh, params, ics)

        if output_write_shp:
            print('Writing shp file...')
            write_shp(base_dir + output_name + '_USM.shp', mesh, params, ics)


        print('Saving mesh to file ' + output_name + '.mesh')
        with open(user_output_dir + output_name + '.mesh', 'w') as outfile:
            dump_json(mesh, outfile)

        # pdb.set_trace()
        if write_param_per_file:
            print('Saving parameters to individual files')
            for key, data in params.items():
                out_path = user_output_dir + f'{output_name}.param.{key}.json'
                with open(out_path, 'w') as outfile:
                    dump_json({key: data}, outfile)
        else:
            print('Saving parameters to file ' + output_name + '.param')
            with open(user_output_dir + output_name + '.param', 'w') as outfile:
                dump_json(params, outfile)

        print('Saving initial conditions  to file ' + output_name + '.ic')
        with open(user_output_dir + output_name + '.ic', 'w') as outfile:
            dump_json(ics, outfile)
        print(f'Export took {time.perf_counter() - export_started:.2f} s')
        print(f'Total run took {time.perf_counter() - run_started:.2f} s')
        print('Done')


def run_parameterization(mesh, parameter_args, base_dir, exec_str, nranks, configfile, output_dir=None):
    """Launch workers once against shared arrays and collect only this run's results."""
    from .mesher_utls.MPI_do_parameterize import rank_face_range

    worker = os.path.join(os.path.dirname(mesher_utls.__file__), 'MPI_do_parameterize.py')
    configfile = os.path.abspath(configfile)
    if nranks < 1:
        raise ValueError('Parameterization requires at least one MPI rank')
    input_started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='mesher-parameters-', dir=os.path.abspath(base_dir)) as work:
        verts_path = os.path.join(work, 'vertices.npy')
        elems_path = os.path.join(work, 'elements.npy')
        def shared_path(key, path, dtype):
            array = mesh['mesh'][key]
            if isinstance(array, np.memmap):
                # Preparation arrays are entire files, not rank slices.
                on_disk = np.load(array.filename, mmap_mode='r')
                if (not isinstance(array.base, np.ndarray) and on_disk.shape == array.shape
                        and on_disk.strides == array.strides and on_disk.dtype == array.dtype):
                    array.flush()
                    return str(array.filename)
            copy_array(path, np.asarray(array, dtype=dtype).reshape(-1, 3))
            return path
        verts_path = shared_path('vertex', verts_path, float)
        elems_path = shared_path('elem', elems_path, np.int64)
        manifest = dict(parameter_args, verts_path=verts_path, elems_path=elems_path,
                        output_dir=work, nranks=nranks)
        if output_dir is not None:
            names = list(dict.fromkeys([*parameter_args['parameter_files'], 'area', 'id']))
            columns = {key: disk_array(os.path.join(output_dir, f'parameter-{i}.npy'),
                                      (len(mesh['mesh']['elem']),), np.int64 if key == 'id' else float)
                       for i, key in enumerate(names)}
            manifest['result_paths'] = {key: str(value.filename) for key, value in columns.items()}
        inputs = os.path.join(work, 'inputs.pickle')
        with open(inputs, 'wb') as stream:
            cloudpickle.dump(manifest, stream)

        inputs_seconds = time.perf_counter() - input_started
        workers_started = time.perf_counter()
        if exec_str is not None:
            command = f'{exec_str} {shlex.join([worker, inputs, "False", configfile])}'
            print(command, flush=True)
            subprocess.check_call(command, shell=True, cwd=os.getcwd())
        else:
            comm = MPI.COMM_SELF.Spawn(sys.executable, args=[worker, inputs, 'True', configfile],
                                       maxprocs=nranks)
            comm.Disconnect()

        workers_seconds = time.perf_counter() - workers_started
        collection_started = time.perf_counter()
        result = []
        nfaces = len(mesh['mesh']['elem'])
        for rank in range(nranks):
            start, stop = rank_face_range(nfaces, rank, nranks)
            with open(os.path.join(work, f'result_{rank}.pickle'), 'rb') as stream:
                batch = cloudpickle.load(stream)
            if output_dir is not None:
                if batch != {'rank': rank, 'start': start, 'stop': stop, 'complete': True}:
                    raise ValueError(f'Invalid parameterization completion for rank {rank}')
                for lo, hi in batches(stop-start):
                    if not np.array_equal(columns['id'][start+lo:start+hi], np.arange(start+lo, start+hi)):
                        raise ValueError(f'Invalid parameterization IDs for rank {rank}')
                continue
            if ((batch['rank'], batch['start'], batch['stop']) != (rank, start, stop)
                    or len(batch['triangles']) != stop-start
                    or any(row['id'] != index for index, row in enumerate(batch['triangles'], start))):
                raise ValueError(f'Invalid parameterization results for rank {rank}')
            result.extend(batch['triangles'])
        print(f'Parameter inputs: {inputs_seconds:.2f} s; '
              f'MPI workers: {workers_seconds:.2f} s; '
              f'result collection: {time.perf_counter() - collection_started:.2f} s', flush=True)
        return columns if output_dir is not None else result


def find_gdal_prefix():
    # we need to make sure we pickup the right paths to all the gdal scripts
    gdal_prefix = ''
    try:
        gdal_prefix = subprocess.run(["gdal-config", "--prefix"], stdout=subprocess.PIPE, env=os.environ, check=True).stdout.decode()
        gdal_prefix = gdal_prefix.replace('\n', '')
        gdal_prefix += '/bin/'
    except:
        raise BaseException(""" ERROR: Could not find gdal-config, please ensure it is installed and on $PATH """)
    return gdal_prefix


def read_config(configfile):
    # Load in configuration file as module
    X = importlib.machinery.SourceFileLoader('config', configfile)
    X = X.load_module()
    # get config file dire
    cwd = os.path.join(os.getcwd(), os.path.dirname(configfile))
    os.chdir(cwd)

    if not hasattr(X, 'dem_filename'):
        raise Exception('dem_filename must be specified in the configuration file')

    if not hasattr(X, 'max_area'):
        raise Exception('max_area must be specified in the configuration file')

    dem_filename = X.dem_filename.strip()
    max_area = X.max_area

    # load any user given parameter files
    parameter_files = {}
    if hasattr(X, 'parameter_files'):
        parameter_files = X.parameter_files

        for key, data in parameter_files.items():
            if isinstance(data['file'], list):
                for i in range(len(data['file'])):
                    parameter_files[key]['file'][i] = parameter_files[key]['file'][i].strip()

    # initial conditions to apply to the triangles with
    initial_conditions = {}
    if hasattr(X, 'initial_conditions'):
        initial_conditions = X.initial_conditions

        for key, data in initial_conditions.items():
            initial_conditions[key]['file'] = initial_conditions[key]['file'].strip()
    # any extra user-specific constraints. E.g., streams
    constraints = {}
    if hasattr(X, 'constraints'):
        constraints = X.constraints

        for key, data in constraints.items():
            constraints[key]['file'] = constraints[key]['file'].strip()
    # simplify applies a simplification routine in GDAL to the outer boundary such
    # that the error between the simplified geom and the original is no more than simplify_tol.
    # This allows the mesh generator more lee way in not creating small triangles allong the boundary
    simplify = False
    simplify_tol = 10
    if hasattr(X, 'simplify'):
        simplify = X.simplify

        if hasattr(X, 'simplify_tol'):
            simplify_tol = X.simplify_tol
    # simplify buffer contracts the outer geometry by bufferDist to help avoid creating small triangles
    bufferDist = -10  # default to -10, will only trigger is simplify is set to T
    if hasattr(X, 'simplify_buffer'):
        bufferDist = X.simplify_buffer
    # Can set simplify_buffer to 0, but we can also just disable this here
    no_simplify_buffer = False
    if hasattr(X, 'no_simplify_buffer'):
        no_simplify_buffer = X.no_simplify_buffer
    if bufferDist > 0:
        print('Buffer must be < 0 as we need to shrink the extent.')
        exit(-1)
    # Enable lloyd iterations "The goal of this mesh optimization is to improve the angles inside the mesh,
    # and make them as close as possible to 60 degrees." 100 iterations is a suggested amount
    # https://doc.cgal.org/latest/Mesh_2/index.html#secMesh_2_optimization
    lloyd_itr = 0
    if hasattr(X, 'lloyd_itr'):
        lloyd_itr = X.lloyd_itr
    # maximum tolerance, as measured by the error metric, between the triangle and the elevation raster.
    # -1 skips the tolerance check -- useful for producing uniform triangles
    max_tolerance = None
    if hasattr(X, 'max_tolerance'):
        max_tolerance = X.max_tolerance
    # Choice of error metric.
    # Can be RMSE or tolerance
    errormetric = 'rmse'
    if hasattr(X, 'errormetric'):
        errormetric = X.errormetric
    # if a mesh was already generated, and only applying a new parametrization is required,
    # enabling this skips the mesh generation step
    reuse_mesh = False
    if hasattr(X, 'reuse_mesh'):
        reuse_mesh = X.reuse_mesh

    # uses GDAL to fill holes
    fill_holes = False
    if hasattr(X, 'fill_holes'):
        fill_holes = X.fill_holes

    mesher_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin', 'meshercxx')

    # look for MESHER_EXE as an environment variable. Defining the mesher path in the config file takes precedenc
    # over this
    using_mesher_environ = False
    try:
        mesher_path = os.environ['MESHER_EXE']
        using_mesher_environ = True
    except KeyError as E:
        pass
    # path to mesher executable
    if hasattr(X, 'mesher_path'):
        mesher_path = X.mesher_path

        if using_mesher_environ:
            warnings.warn(
                "Warning: mesher binary path defined in env var and in configuration file. Using the mesher path from "
                "the configuration file")

    # enable verbose output for debugging
    verbose = False
    if hasattr(X, 'verbose'):
        verbose = X.verbose
    mesher_debug = False
    if hasattr(X, 'mesher_debug'):
        mesher_debug = X.mesher_debug
    # Default to exactextract for speed; allow users to opt out to retain
    # the slower per-pixel parameterization semantics.
    use_exactextract = True
    if hasattr(X, 'use_exactextract'):
        use_exactextract = X.use_exactextract
    write_param_per_file = False
    if hasattr(X, 'write_param_per_file'):
        write_param_per_file = X.write_param_per_file
    user_output_dir = cwd + os.path.sep
    # output to the specific directory, instead of the root dir of the calling python script
    if hasattr(X, 'user_output_dir'):
        user_output_dir += X.user_output_dir
    else:
        # use the config filename as output path
        (file, ext) = os.path.splitext(os.path.basename(configfile))
        user_output_dir += file
    if user_output_dir[-1] is not os.path.sep:
        user_output_dir += os.path.sep
    # should we write a shape file of the USM, pretty costly on the big meshes
    output_write_shp = True
    if hasattr(X, 'write_shp'):
        output_write_shp = X.write_shp
    output_write_vtu = True
    if hasattr(X, 'write_vtu'):
        output_write_vtu = X.write_vtu
    # Use the input file's projection.
    # This is useful for preserving a UTM input. Does not work if the input file is geographic.
    use_input_prj = True
    if hasattr(X, 'use_input_prj'):
        use_input_prj = X.use_input_prj
    # Do smoothing of the input DEM to create a more smooth mesh. This can help if the DEM quality is poor or if
    # triangles close to the elevation raster cell size is required
    do_smoothing = False
    if hasattr(X, 'do_smoothing'):
        do_smoothing = X.do_smoothing
    # Smoothing factor for above option.
    scaling_factor = 2.0
    if hasattr(X, 'smoothing_scaling_factor'):
        scaling_factor = X.smoothing_scaling_factor
    # number of iterations to smooth over, each smoothing using cubic spline,
    # and resamples by iter * scaling_factor for each iteration
    max_smooth_iter = 1
    if hasattr(X, 'max_smooth_iter'):
        max_smooth_iter = X.max_smooth_iter
    # use the convex combination of weights method
    user_no_weights = False
    weight_threshold = 0.5
    if hasattr(X, 'weight_threshold'):
        weight_threshold = X.weight_threshold
    if hasattr(X, 'use_weights'):
        if not X.use_weights:
            user_no_weights = True
    wkt_out = "PROJCS[\"North_America_Albers_Equal_Area_Conic\"," \
              "     GEOGCS[\"GCS_North_American_1983\"," \
              "         DATUM[\"North_American_Datum_1983\"," \
              "             SPHEROID[\"GRS_1980\",6378137,298.257222101]]," \
              "         PRIMEM[\"Greenwich\",0]," \
              "         UNIT[\"Degree\",0.017453292519943295]]," \
              "     PROJECTION[\"Albers_Conic_Equal_Area\"]," \
              "     PARAMETER[\"False_Easting\",0]," \
              "     PARAMETER[\"False_Northing\",0]," \
              "     PARAMETER[\"longitude_of_center\",-96]," \
              "     PARAMETER[\"Standard_Parallel_1\",20]," \
              "     PARAMETER[\"Standard_Parallel_2\",60]," \
              "     PARAMETER[\"latitude_of_center\",40]," \
              "     UNIT[\"Meter\",1]," \
              "     AUTHORITY[\"EPSG\",\"102008\"]]"
    if hasattr(X, 'wkt_out'):
        wkt_out = X.wkt_out
    # use a custom extent instead of the entire input DEM to define the meshing region
    # extent = [xmin, ymin, xmax, ymax] in source SRS
    extent = None
    if hasattr(X, 'extent'):
        extent = X.extent
    # Clip to a shape file
    clip_to_shp = None
    if hasattr(X, 'clip_to_shp'):
        clip_to_shp = X.clip_to_shp

        if not os.path.exists(clip_to_shp):
            raise Exception(f'Clipping shape file is not valid. Path given was\n {clip_to_shp}')
    if clip_to_shp and extent:
        raise Exception('Cannot specify both extent and a shape file to clip to')
    nworkers = os.cpu_count() or 1
    # on linux we can ensure that we respect cpu affinity
    if 'sched_getaffinity' in dir(os):
        nworkers = len(os.sched_getaffinity(0))
    if hasattr(X, 'nworkers'):
        nworkers = X.nworkers
    # GDAL workers are for the initial regularize call where multiple gdalwarp processes are started
    # these are memory heavy so this can be limited if required
    nworkers_gdal = nworkers
    if hasattr(X, 'nworkers_gdal'):
        nworkers_gdal = X.nworkers_gdal
    try:
        nworkers_gdal = os.environ['MESHER_NWORKERS_GDAL']

        if hasattr(X, 'nworkers_gdal'):
            print('Warning: Overridding configfile nworkers_gdal with environment variable MESHER_NWORKERS_GDAL')
    except KeyError as E:
        pass
    try:
        nworkers = os.environ['MESHER_NWORKERS']

        if hasattr(X, 'nworkers'):
            print('Warning: Overridding configfile nworkers with environment variable MESHER_NWORKERS')
    except KeyError as E:
        pass
    nworkers = int(nworkers)
    nworkers_gdal = int(nworkers_gdal)
    try:
        cache = os.environ['GDAL_CACHEMAX']
        print(
            'Warning: GDAL_CACHEMAX={0}. Ensure this is reasonable for the number of parallel processes specified'.format(
                cache))
    except:
        # the user doesn't have the GDAL_CACHEMAX defined, which is good. Set it to a reasonable value
        os.environ['GDAL_CACHEMAX'] = "{0}%".format(30.0 / nworkers)
    print('Using {0} CPUs and per-CPU memory cache = {1}'.format(nworkers, os.environ['GDAL_CACHEMAX']))

    MPI_exec_str = None
    if hasattr(X, 'MPI_exec_str'):
        MPI_exec_str = X.MPI_exec_str

    # on macos M1, the 2 efficiency cores don't seem to be targetable by OpenMPI?
    MPI_nworkers = nworkers
    if hasattr(X, 'MPI_nworkers'):
        MPI_nworkers = X.MPI_nworkers
    elif not hasattr(X, 'MPI_nworkers') and MPI_exec_str:
        raise RuntimeError('If MPI_exec_str is provided, then MPI_nworkers must also be provided')

    _positive_integer(MPI_nworkers, 'MPI_nworkers')
    mpi_vendor = MPI.get_vendor()[0]
    print(f'Detected MPI = {mpi_vendor}')

    # If we are using OpenMPI we can use a blocking spawn, otherwise, we need to use mpirun
    if MPI_exec_str is None and mpi_vendor != 'Open MPI':
        MPI_exec_str = f'mpirun -n {MPI_nworkers} {sys.executable}'

    mpi_mesh = False
    if hasattr(X, 'mpi_mesh'):
        mpi_mesh = X.mpi_mesh

    return X, bufferDist, clip_to_shp, constraints, dem_filename, do_smoothing, errormetric, extent, fill_holes, \
        initial_conditions, lloyd_itr, max_area, max_smooth_iter, max_tolerance, mesher_path, no_simplify_buffer, \
        nworkers, nworkers_gdal, output_write_shp, output_write_vtu, parameter_files, reuse_mesh, scaling_factor, \
        simplify, simplify_tol, use_input_prj, user_no_weights, user_output_dir, verbose, weight_threshold, \
        wkt_out, MPI_exec_str, MPI_nworkers, mpi_mesh, \
        mesher_debug, use_exactextract, write_param_per_file


def regularize_inputs(parameter_dir, gdal_prefix, input_files, pixel_height, pixel_width, srs_out,
                      xmax, xmin, ymax, ymin, fill_holes, MPI_exec_str, MPI_nworkers, reuse_projected_parameters):
    # ensure all the weights sum to 1
    total_weights = 0
    use_weights = False
    param_args = []
    if not input_files:
        return total_weights, use_weights
    os.makedirs(parameter_dir, exist_ok=True)
    for key, data in input_files.items():
        if 'weight' in data:
            total_weights += data['weight']
            use_weights = True
        if 'classifier' in data:
            input_files[key]['classifier'] = cloudpickle.dumps(input_files[key]['classifier'])

        param_args.append((parameter_dir, data, gdal_prefix, key, pixel_height, pixel_width,
                           srs_out.ExportToProj4(), xmax, xmin, ymax, ymin, fill_holes, reuse_projected_parameters))

    with open('pickled_param_args.pickle', 'wb') as f:
        cloudpickle.dump(param_args, f)

    MPI_regularize_path = os.path.join(os.path.dirname(mesher_utls.__file__),
                                                      'MPI_regularize_inputs.py')
    if MPI_exec_str is not None:
        exec_str = f"""{MPI_exec_str} {MPI_regularize_path} pickled_param_args.pickle False {MPI_nworkers}"""
        print(exec_str)
        subprocess.check_call([exec_str], shell=True, cwd=os.getcwd())
    else:
        comm = MPI.COMM_SELF.Spawn(sys.executable,
                                   args=[MPI_regularize_path,
                                         'pickled_param_args.pickle', 'True'],
                                   maxprocs=MPI_nworkers)
        comm.Disconnect()
    os.remove('pickled_param_args.pickle')

    ret = []
    for file in glob.glob('pickled_param_args_rets_*.pickle'):
        with open(file, 'rb') as f:
            ret.append(cloudpickle.load(f))
        os.remove(file)

    for w in ret:
        if len(w) == 0:
            continue

        # it comes as a [ {} ]

        # r will a list of dicts, len(r) is the amount of work that MPI rank did
        for r in w:
            key = r['key']
            input_files[key]['filename'] = r['filename']
            input_files[key]['file'] = None

            if not isinstance(input_files[key]['method'], list):
                input_files[key]['method'] = [input_files[key]['method']]

    return total_weights, use_weights


def extract_point(raster, mx, my):
    # Convert from map to pixel coordinates.
    # Only works for geotransforms with no rotation.

    rb = raster.GetRasterBand(1)
    gt = raster.GetGeoTransform()

    px = int((mx - gt[0]) / gt[1])  # x pixel
    py = int((my - gt[3]) / gt[5])  # y pixel

    # boundary verticies from Triangle can end up outside of the domain by 1 pixel.
    # if we adjusted back by  1 pixel to get the dz value, it's no problem and still gives a good boundary
    if px == raster.RasterXSize:
        px = px - 1
    if py == raster.RasterYSize:
        py = py - 1
    try:
        mz = rb.ReadAsArray(px, py, 1, 1)
    except RuntimeError:
        mz = rb.GetNoDataValue()

    if mz is None:
        return rb.GetNoDataValue()
    mz = float(np.asarray(mz).item())

    if mz == rb.GetNoDataValue():

        # look at the surrounding 8 cells
        #  (px-1, py - 1)  (px, py - 1)  (px + 1,py - 1)
        #  (px-1, py    )    (px, py)      (px + 1, py   )
        #  (px-1, py + 1)  (px, py + 1)  (px + 1, py + 1)

        coords = [
            (px - 1, py - 1), (px, py - 1), (px + 1, py - 1),
            (px - 1, py),     (px, py    ), (px + 1, py    ),
            (px - 1, py + 1), (px, py + 1), (px + 1, py + 1)
        ]

        zs = []
        for c in coords:
            try:
                z1 = rb.ReadAsArray(c[0], c[1], 1, 1)
                zs.append(z1)
            except RuntimeError as e:
                pass  # out of bounds, pass

        z = [x for x in zs if x != rb.GetNoDataValue() and x is not None ]

        if len(z) == 0:
            # print 'Warning: The point (%s,%s) and its 8-neighbours lies outside of the DEM domain' % (mx, my)
            return rb.GetNoDataValue()

        mz = float(np.mean(z))

    return mz


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return value


def compute_tile_grid(nranks, width, height, *, num_tiles=None):
    """Balance worker utilization with the physical perimeter of each tile.

    A square tile has shape cost 1. Multiply that cost by ranks/tiles so a
    slightly rectangular grid can use more workers without choosing long thin
    strips just to fill a prime-sized allocation. This is a geometry heuristic,
    not an estimate of raster-dependent refinement work. An explicit num_tiles
    chooses the best-shaped exact grid independently of the worker count.
    """
    _positive_integer(nranks, 'MPI_nworkers')
    if not all(np.isfinite(v) and v > 0 for v in (width, height)):
        raise ValueError('MPI domain width and height must be positive and finite')
    if num_tiles is not None:
        _positive_integer(num_tiles, 'mpi_num_tiles')
        # An explicit count must be exact. Of its factor pairs, find the
        # grid with the most nearly square physical tiles.
        candidates = ((short, num_tiles // short) if width >= height else (num_tiles // short, short)
                      for short in range(1, math.isqrt(num_tiles) + 1) if num_tiles % short == 0)
        return min(candidates, key=lambda grid: abs(math.log(width / grid[1]) - math.log(height / grid[0])))

    best = None
    for rows in range(1, nranks + 1):
        for cols in range(1, nranks // rows + 1):
            count = rows * cols
            aspect = (width / cols) / (height / rows)
            shape_cost = (math.sqrt(aspect) + 1 / math.sqrt(aspect)) / 2
            score = shape_cost * nranks / count
            candidate = (score, -count, abs(rows-cols), rows, cols)
            if best is None or candidate < best:
                best = candidate
    _, _, _, rows, cols = best
    return rows, cols


def mpi_final_npz_name(base_name, rows, cols):
    # Pairwise stages halve columns and rows in turn, carrying odd groups.
    stages = (cols - 1).bit_length() + (rows - 1).bit_length()
    suffix = f'stage{stages-1}_r0_c0' if stages else 'tile_0'
    return f'{base_name}_mpi_{suffix}.npz'


def build_mesh_from_arrays(verts_xy, tris, src_ds, dem, srs, is_geographic, verbose, directory=None):
    """Keep one numeric representation; sample and repair elevations in place."""
    verts_xy = np.asarray(verts_xy, dtype=float)
    tris = integer_rows(tris, 'elem')
    def allocate(name, shape, dtype):
        return (np.empty(shape, dtype=dtype) if directory is None else
                disk_array(os.path.join(directory, name + '.npy'), shape, dtype))
    elements = allocate('elements', tris.shape, np.int64)
    for lo, hi in batches(len(tris)):
        elements[lo:hi] = tris[lo:hi]
    neighbours = build_neighbours(elements, len(verts_xy), verts_xy, chunk_size=32768,
        output_path=None if directory is None else os.path.join(directory, 'neighbours.npy'),
        temporary_directory=directory)
    with tempfile.TemporaryDirectory(dir=directory) as work:
        check_boundary_geometry(verts_xy, elements, neighbours, work)
    vertices = allocate('vertices', (len(verts_xy), 3), float)
    invalid = allocate('invalid-elevations', (len(verts_xy),), bool)
    for lo, hi in batches(len(vertices)):
        vertices[lo:hi, :2] = verts_xy[lo:hi, :2]
        for i in range(lo, hi):
            z = extract_point(src_ds, float(vertices[i, 0]), float(vertices[i, 1]))
            invalid[i] = z is None or z == dem.GetNoDataValue() or not np.isfinite(z)
            vertices[i, 2] = np.nan if invalid[i] else z
    # Preserve the original ordered, one-pass propagation through face vertices.
    for lo, hi in batches(len(elements)):
        for tri in elements[lo:hi]:
            for vid in tri:
                if invalid[vid]:
                    valid = tri[~invalid[tri]]
                    if len(valid):
                        vertices[vid, 2] = vertices[valid, 2].mean()
                        invalid[vid] = False
    count = sum(np.count_nonzero(invalid[lo:hi]) for lo, hi in batches(len(invalid)))
    if count:
        raise RuntimeError(f'Length of invalid nodes after correction= {count}. '
                           'An entire triangle may lie outside the DEM; try reducing simplify_tol.')
    m = dict(vertex=vertices, elem=elements, neigh=neighbours, nvertex=len(vertices),
             nelem=len(elements), is_geographic=is_geographic)
    if not is_geographic:
        m.update(proj4=srs.ExportToProj4(), UTM_zone=srs.GetUTMZone())
    return {'mesh': m}


def load_mesh_from_mesher_files(base_dir, base_name, src_ds, dem, srs, is_geographic, verbose, directory=None):
    """Read Triangle text files directly into numeric backing arrays."""
    def rows(path):
        with open(path) as stream:
            for line in stream:
                fields = line.split('#', 1)[0].split()
                if fields:
                    yield fields
    with tempfile.TemporaryDirectory(prefix='serial-arrays-', dir=directory) as work:
        nodes = rows(base_dir + 'PLGS' + base_name + '.1.node')
        nv = int(next(nodes)[0])
        vertices = disk_array(os.path.join(work, 'vertices.npy'), (nv, 3))
        count = 0
        for count, fields in enumerate(nodes, 1):
            if int(fields[0]) != count or count > nv:
                raise MeshTopologyError('Unexpected node ID in Triangle output')
            vertices[count-1] = (float(fields[1]), float(fields[2]), 0)
        if count != nv:
            raise MeshTopologyError('Truncated Triangle node file')
        faces = rows(base_dir + 'PLGS' + base_name + '.1.ele')
        nt = int(next(faces)[0])
        elements = disk_array(os.path.join(work, 'elements.npy'), (nt, 3), np.int64)
        count = 0
        for count, fields in enumerate(faces, 1):
            if int(fields[0]) != count or count > nt:
                raise MeshTopologyError('Unexpected face ID in Triangle output')
            elements[count-1] = [int(v)-1 for v in fields[1:4]]
        if count != nt:
            raise MeshTopologyError('Truncated Triangle element file')
        return build_mesh_from_arrays(vertices, elements, src_ds, dem, srs, is_geographic,
                                      verbose, directory=directory)


def _run_mpi_mesh_tasks(tasks, MPI_exec_str, MPI_nworkers, worker, label):
    if not tasks:
        return
    debug = any(task.get('mesher_debug', False) for task in tasks)
    worker_path = os.path.join(os.path.dirname(mesher_utls.__file__), worker)
    with StageProgress(tasks, label, debug=debug) as progress:
        arguments = str(progress.directory / 'tasks.pickle')
        with open(arguments, 'wb') as stream:
            cloudpickle.dump(progress.tasks, stream)
        if MPI_exec_str is not None:
            exec_str = f'{MPI_exec_str} {shlex.quote(worker_path)} {shlex.quote(arguments)} False {MPI_nworkers}'
            if debug:
                print(exec_str, flush=True)
            subprocess.check_call([exec_str], shell=True, cwd=os.getcwd())
        else:
            comm = MPI.COMM_SELF.Spawn(sys.executable,
                                       args=[worker_path, arguments, 'True'],
                                       maxprocs=min(MPI_nworkers, len(tasks)))
            comm.Disconnect()


def run_mpi_merge_stage(tasks, MPI_exec_str, MPI_nworkers, stage_label=None, total_merges=None):
    if not tasks:
        return
    if stage_label is not None and any(task.get('mesher_debug', False) for task in tasks):
        remaining = total_merges - len(tasks) if total_merges is not None else None
        if remaining is not None:
            print(f'MPI merge {stage_label}: {len(tasks)} merges (remaining after stage: {remaining})')
        else:
            print(f'MPI merge {stage_label}: {len(tasks)} merges')

    _run_mpi_mesh_tasks(tasks, MPI_exec_str, MPI_nworkers, 'MPI_merge_tiles.py',
                        f'MPI merge {stage_label}' if stage_label else 'MPI merge')


def mpi_remesh_options_from_config(config):
    """Validate MPI meshing settings before preparing or replacing any inputs."""
    if getattr(config, 'mpi_merge_method', 'remesh') != 'remesh':
        raise ValueError('Only mpi_merge_method="remesh" is supported; shared-edge meshing has been retired')
    for name in ('mpi_merge_snap_tol', 'mpi_shared_edge_spacing', 'mpi_shared_edge_constraints'):
        if getattr(config, name, None) is not None:
            raise ValueError(f'{name} belongs to the retired shared-edge implementation; remove it from the configuration')
    if getattr(config, 'mpi_global_lloyd', 0) != 0:
        raise ValueError('Distributed MPI smoothing has been retired; remove mpi_global_lloyd or set it to 0')
    if getattr(config, 'mpi_lloyd_boundary_tol', None) is not None:
        raise ValueError('mpi_lloyd_boundary_tol belongs to the retired distributed smoother; remove it from the configuration')
    if getattr(config, 'lloyd_itr', 0) != 0:
        raise ValueError('lloyd_itr must be 0 with MPI remeshing; smoothing does not preserve the validated patch criteria')
    return {'merge_method': 'remesh',
            'remesh_band_width': getattr(config, 'mpi_remesh_band_width', None),
            'remesh_max_attempts': getattr(config, 'mpi_remesh_max_attempts', 12),
            'remesh_max_vertices': getattr(config, 'mpi_remesh_max_vertices', 500000),
            'remesh_max_active_faces': getattr(config, 'mpi_remesh_max_active_faces', 2000000)}


def validate_mpi_mesh_metadata(npz_path):
    """Only completed cavity meshes can be reused without remeshing geometry."""
    metadata_path = os.path.splitext(npz_path)[0] + '.json'
    try:
        with open(metadata_path) as stream:
            metadata = json.load(stream)
    except (OSError, ValueError) as exc:
        raise ValueError(f'Cannot verify MPI mesh metadata at {metadata_path}; regenerate with reuse_mesh=False') from exc
    if not isinstance(metadata, dict) or metadata.get('merge_method') != 'remesh':
        raise ValueError('Cannot reuse an old shared-edge MPI mesh; regenerate with reuse_mesh=False')


def run_mpi_meshing(base_dir, base_name, xmin, ymin, xmax, ymax, mesher_path,
                    outer_polygon_shp, constraints, parameter_files, initial_conditions, max_area, min_area,
                    max_tolerance, errormetric, use_weights, topo_weight,
                    weight_threshold, is_geographic, MPI_exec_str, MPI_nworkers,
                    mesher_debug=False, remesh_options=None, *, num_tiles=None):
    remesh_options = dict(remesh_options or {})
    if remesh_options.setdefault('merge_method', 'remesh') != 'remesh':
        raise ValueError('Only mpi_merge_method="remesh" is supported; shared-edge meshing has been retired')
    if is_geographic:
        raise ValueError('MPI cavity remeshing requires a projected working CRS')
    for key, minimum, default in (('remesh_max_attempts', 1, 12), ('remesh_max_vertices', 3, 500000),
                                  ('remesh_max_active_faces', 1, 2000000)):
        value = remesh_options.get(key, default)
        if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
            raise ValueError(f'mpi_{key} must be an integer >= {minimum}')
    rows, cols = compute_tile_grid(MPI_nworkers, xmax-xmin, ymax-ymin, num_tiles=num_tiles)
    layout_rows, layout_cols = rows, cols

    tile_width = (xmax - xmin) / cols
    tile_height = (ymax - ymin) / rows
    min_edge = min(tile_width, tile_height)
    band_width = remesh_options.get('remesh_band_width')
    if band_width is None:
        band_width = min(2 * math.sqrt(max_area), .2 * min_edge)
    if not np.isfinite(band_width) or band_width <= 0:
        raise ValueError('mpi_remesh_band_width must be positive and finite')
    remesh_options['remesh_band_width'] = float(band_width)
    min_edge_threshold = 1.25 * band_width
    if min_edge <= min_edge_threshold:
        msg = (
            f'The {rows}x{cols} MPI tile grid has tiles '
            f'with min edge {min_edge:.3f} <= threshold {min_edge_threshold:.3f}. '
            'The remesh band is too wide relative to the tiles. '
            'Use fewer tiles or reduce mpi_remesh_band_width.'
        )
        raise RuntimeError(msg)

    tile_count = rows * cols
    print(f'MPI tile grid: {rows} rows x {cols} columns; '
          f'{tile_count} tiles, up to {min(tile_count, MPI_nworkers)} active tile workers; '
          f'tile dimensions={tile_width:.6g} x {tile_height:.6g}', flush=True)


    constraint_files = [v['filename'] for v in constraints.values() if 'filename' in v]

    tile_args = []
    # Shared coordinates are computed once, preventing roundoff differences in
    # neighbouring core boxes and making ownership an exact geometric contract.
    x_edges = np.linspace(xmin, xmax, cols + 1)
    y_edges = np.linspace(ymin, ymax, rows + 1)
    for tile_id in range(tile_count):
        row, col = divmod(tile_id, cols)

        tile_xmin, tile_xmax = x_edges[col:col + 2]
        tile_ymin, tile_ymax = y_edges[rows - row - 1:rows - row + 1]

        tile_bbox = [tile_xmin, tile_ymin, tile_xmax, tile_ymax]
        tile_prefix = f"{base_name}_mpi_tile_{tile_id}"

        tile_args.append({
            'tile_id': tile_id,
            'tile_bbox': tile_bbox,
            'band_width': band_width,
            'base_dir': base_dir,
            'base_name': base_name,
            'tile_prefix': tile_prefix,
            'mesher_path': mesher_path,
            'constraints': constraint_files,
            'parameter_files': parameter_files,
            'initial_conditions': initial_conditions,
            'max_area': max_area,
            'min_area': min_area,
            'max_tolerance': max_tolerance,
            'errormetric': errormetric,
            'use_weights': use_weights,
            'topo_weight': topo_weight,
            'weight_threshold': weight_threshold,
            'is_geographic': is_geographic,
            'dem_path': base_dir + base_name + '_projected.tif',
            'outer_polygon_shp': outer_polygon_shp,
            'mesher_debug': mesher_debug,
            **remesh_options
        })

    _run_mpi_mesh_tasks(tile_args, MPI_exec_str, MPI_nworkers, 'MPI_mesh_tiles.py', 'MPI tiles')

    tile_grid = []
    for row in range(rows):
        row_tiles = []
        for col in range(cols):
            tile_id = row * cols + col
            tile_prefix = f"{base_name}_mpi_tile_{tile_id}"
            row_tiles.append({
                'npz': base_dir + tile_prefix + '.npz',
                'meta': base_dir + tile_prefix + '.json'
            })
        tile_grid.append(row_tiles)

    stage = 0
    total_merges = 0
    tmp_rows = rows
    tmp_cols = cols
    while tmp_rows > 1 or tmp_cols > 1:
        if tmp_cols > 1:
            total_merges += tmp_rows * (tmp_cols // 2)
            tmp_cols = (tmp_cols + 1) // 2
        if tmp_rows > 1:
            total_merges += (tmp_rows // 2) * tmp_cols
            tmp_rows = (tmp_rows + 1) // 2

    completed_merges = 0
    while rows > 1 or cols > 1:
        if cols > 1:
            tasks = []
            new_grid = []
            for row in range(rows):
                new_row = []
                col = 0
                while col < cols:
                    if col + 1 < cols:
                        left = tile_grid[row][col]
                        right = tile_grid[row][col + 1]
                        out_prefix = f"{base_name}_mpi_stage{stage}_r{row}_c{col // 2}"
                        tasks.append({
                            'tile_a': left,
                            'tile_b': right,
                            'out_prefix': base_dir + out_prefix,
                            'mesher_path': mesher_path,
                            'constraints': constraint_files,
                            'parameter_files': parameter_files,
                            'initial_conditions': initial_conditions,
                            'max_area': max_area,
                            'min_area': min_area,
                            'max_tolerance': max_tolerance,
                            'errormetric': errormetric,
                            'use_weights': use_weights,
                            'topo_weight': topo_weight,
                            'weight_threshold': weight_threshold,
                            'is_geographic': is_geographic,
                            'dem_path': base_dir + base_name + '_projected.tif',
                            'outer_polygon_shp': outer_polygon_shp,
                            'mesher_debug': mesher_debug,
                            **remesh_options
                        })
                        new_row.append({
                            'npz': base_dir + out_prefix + '.npz',
                            'meta': base_dir + out_prefix + '.json'
                        })
                        col += 2
                    else:
                        new_row.append(tile_grid[row][col])
                        col += 1
                new_grid.append(new_row)

            run_mpi_merge_stage(
                tasks, MPI_exec_str, MPI_nworkers,
                stage_label=f'stage{stage} cols',
                total_merges=total_merges - completed_merges
            )
            completed_merges += len(tasks)
            tile_grid = new_grid
            cols = len(tile_grid[0])
            stage += 1

        if rows > 1:
            tasks = []
            new_grid = []
            row = 0
            new_row_idx = 0
            while row < rows:
                if row + 1 < rows:
                    merged_row = []
                    for col in range(cols):
                        top = tile_grid[row][col]
                        bottom = tile_grid[row + 1][col]
                        out_prefix = f"{base_name}_mpi_stage{stage}_r{new_row_idx}_c{col}"
                        tasks.append({
                            'tile_a': top,
                            'tile_b': bottom,
                            'out_prefix': base_dir + out_prefix,
                            'mesher_path': mesher_path,
                            'constraints': constraint_files,
                            'parameter_files': parameter_files,
                            'initial_conditions': initial_conditions,
                            'max_area': max_area,
                            'min_area': min_area,
                            'max_tolerance': max_tolerance,
                            'errormetric': errormetric,
                            'use_weights': use_weights,
                            'topo_weight': topo_weight,
                            'weight_threshold': weight_threshold,
                            'is_geographic': is_geographic,
                            'dem_path': base_dir + base_name + '_projected.tif',
                            'outer_polygon_shp': outer_polygon_shp,
                            'mesher_debug': mesher_debug,
                            **remesh_options
                        })
                        merged_row.append({
                            'npz': base_dir + out_prefix + '.npz',
                            'meta': base_dir + out_prefix + '.json'
                        })
                    new_grid.append(merged_row)
                    row += 2
                    new_row_idx += 1
                else:
                    new_grid.append(tile_grid[row])
                    row += 1
                    new_row_idx += 1

            run_mpi_merge_stage(
                tasks, MPI_exec_str, MPI_nworkers,
                stage_label=f'stage{stage} rows',
                total_merges=total_merges - completed_merges
            )
            completed_merges += len(tasks)
            tile_grid = new_grid
            rows = len(tile_grid)
            stage += 1

    final_path = tile_grid[0][0]['npz']
    layout_path = os.path.join(base_dir, base_name + '_mpi_layout.json')
    temporary = layout_path + '.tmp'
    with open(temporary, 'w') as stream:
        json.dump({'version': 1, 'rows': layout_rows, 'cols': layout_cols,
                   'mpi_nworkers': MPI_nworkers, 'num_tiles': tile_count, 'bounds': [xmin, ymin, xmax, ymax],
                   'final_npz': os.path.basename(final_path)}, stream)
    os.replace(temporary, layout_path)
    return final_path


def compute_mpi_final_npz_path(base_dir, base_name, mpi_nworkers):
    """Locate the completed layout, with compatibility for pre-manifest runs."""
    layout_path = os.path.join(base_dir, base_name + '_mpi_layout.json')
    if os.path.exists(layout_path):
        try:
            with open(layout_path) as stream:
                layout = json.load(stream)
            rows, cols = layout['rows'], layout['cols']
            if (layout.get('version') != 1 or any(type(v) is not int or v < 1 for v in (rows, cols))
                    or layout['final_npz'] != mpi_final_npz_name(base_name, rows, cols)):
                raise ValueError('Invalid layout fields')
        except (OSError, ValueError, TypeError, KeyError) as exc:
            raise ValueError(f'Cannot read completed MPI layout {layout_path}; '
                             'regenerate with reuse_mesh=False') from exc
        return os.path.join(base_dir, layout['final_npz'])
    # Older versions always selected the largest square grid <= the rank count.
    # Use their naming rule only for old cached outputs, never for fresh meshing.
    if not isinstance(mpi_nworkers, int) or isinstance(mpi_nworkers, bool) or mpi_nworkers < 1:
        raise ValueError('MPI_nworkers must be a positive integer')
    side = math.isqrt(mpi_nworkers)
    return os.path.join(base_dir, mpi_final_npz_name(base_name, side, side))


# To remove the dependency on gdal python bindings and thus the scripts, this is a slimmed down version of the
# gdal_polygonize.py code from
# https://github.com/OSGeo/gdal/blob/release/2.4/gdal/swig/python/scripts/gdal_polygonize.py
def _geom_point_count(geom):
    # Total number of coordinate points across all rings of a (multi)polygon or linestring.
    # Useful for logging how much a simplify/buffer step reduced the boundary complexity.
    if geom is None:
        return -1
    count = 0
    stack = [geom]
    while stack:
        g = stack.pop()
        t = g.GetGeometryType()
        if t in (ogr.wkbPolygon, ogr.wkbPolygon25D):
            for i in range(g.GetGeometryCount()):
                stack.append(g.GetGeometryRef(i))
        elif t in (ogr.wkbLinearRing, ogr.wkbLineString, ogr.wkbLineString25D):
            count += g.GetPointCount()
        else:
            for i in range(g.GetGeometryCount()):
                stack.append(g.GetGeometryRef(i))
    return count


def gdal_polygonize(src_filename, mask, dst_filename):
    options = []
    src_band_n = 1
    src_ds = gdal.Open(src_filename)

    srcband = src_ds.GetRasterBand(src_band_n)

    mask_ds = gdal.Open(mask)
    maskband = mask_ds.GetRasterBand(1)

    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(dst_filename)

    srs = None
    if src_ds.GetProjectionRef() != '':
        srs = osr.SpatialReference()

        srs.ImportFromWkt(src_ds.GetProjectionRef())

    dst_layer = dst_ds.CreateLayer('out', geom_type=ogr.wkbPolygon, srs=srs)

    fd = ogr.FieldDefn('DN', ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = 0

    result = gdal.Polygonize(srcband, maskband, dst_layer, dst_field, options)


def write_shp(fname, mesh, parameter_files, initial_conditions):
    # Create the shape file to hold the triangulation

    driver = ogr.GetDriverByName("ESRI Shapefile")

    try:
        os.remove(fname)  # remove if existing
    except OSError:
        pass

    output_usm = driver.CreateDataSource(fname)

    srs_out = osr.SpatialReference()

    srs_out.ImportFromProj4(mesh['mesh']['proj4'])

    layer = output_usm.CreateLayer('mesh', srs_out, ogr.wkbPolygon)

    for key, value in parameter_files.items():
        layer.CreateField(ogr.FieldDefn(key, ogr.OFTReal))

    for key, value in initial_conditions.items():
        layer.CreateField(ogr.FieldDefn(key, ogr.OFTReal))

    for elem in range(mesh['mesh']['nelem']):
        v0 = mesh['mesh']['elem'][elem][0]
        v1 = mesh['mesh']['elem'][elem][1]
        v2 = mesh['mesh']['elem'][elem][2]

        # we need this to do the area calculation
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(mesh['mesh']['vertex'][v0][0], mesh['mesh']['vertex'][v0][1])
        ring.AddPoint(mesh['mesh']['vertex'][v1][0], mesh['mesh']['vertex'][v1][1])
        ring.AddPoint(mesh['mesh']['vertex'][v2][0], mesh['mesh']['vertex'][v2][1])
        ring.AddPoint(mesh['mesh']['vertex'][v0][0],
                      mesh['mesh']['vertex'][v0][1])  # add again to complete the ring.

        # need this for the area calculation
        tpoly = ogr.Geometry(ogr.wkbPolygon)
        tpoly.AddGeometry(ring)

        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(tpoly)

        # for this section,
        # key[0:10] -> if the name is longer, it'll have been truncated when we made the field
        # pdb.set_trace()
        for key, data in parameter_files.items():
            output = data[elem]
            if output is None:
                output = float('nan')
            elif isinstance(output, (list, tuple, np.ndarray)):
                output = output[0] if len(output) else float('nan')
            feature.SetField(key[0:10], float(output))

        for key, data in initial_conditions.items():
            output = data[elem]
            if output is None:
                output = float('nan')
            elif isinstance(output, (list, tuple, np.ndarray)):
                output = output[0] if len(output) else float('nan')
            feature.SetField(key[0:10], float(output))

        layer.CreateFeature(feature)
    output_usm.FlushCache()
    output_usm = None  # close file



if __name__ == "__main__":
    main()
