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
import numpy as np
from functools import partial
import sys
import shutil
import importlib
import vtk
import warnings
import uuid
import cloudpickle
from concurrent import futures
import time

# This reverts to Python <= 3.7 behaviours. It "worked" before...
# Long term it's almost certainly a Bad Idea and should be resolved.
# > Changed in version 3.8: On macOS, the spawn start method is now the default.
# >The fork start method should be considered unsafe as it can lead to crashes of the subprocess. See bpo-33725.
# Setting this to 'spawn' results in the imported config file not being available for the child processes
# Ref https://bugs.python.org/issue33725
import multiprocessing as mp
mp.set_start_method("fork")


gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()  # Enable exception support
osr.UseExceptions()  # Enable exception support


def main():
    # load user configurable parameters here
    # Check user defined configuration file

    if len(sys.argv) == 1:
        print('ERROR: mesher.py requires one argument [configuration file] (i.e. mesher.py Bow.py)')
        return

    # Get name of configuration file/module
    configfile = sys.argv[-1]

    # Load in configuration file as module
    X = importlib.machinery.SourceFileLoader('config', configfile)
    X = X.load_module()

    # get config file dire
    cwd = os.path.join(os.getcwd(), os.path.dirname(configfile))
    os.chdir(cwd)

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
        bufferDist = X.no_simplify_buffer

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

    mesher_path = os.path.dirname(os.path.abspath(__file__)) + '/mesher'

    # uses GDAL to fill holes
    fill_holes = False
    if hasattr(X, 'fill_holes'):
        fill_holes = X.fill_holes


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

    #Clip to a shape file
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
    ########################################################

    # we need to make sure we pickup the right paths to all the gdal scripts
    gdal_prefix = ''
    try:
        gdal_prefix = subprocess.run(["gdal-config", "--prefix"], stdout=subprocess.PIPE).stdout.decode()
        gdal_prefix = gdal_prefix.replace('\n', '')
        gdal_prefix += '/bin/'
    except:
        raise BaseException(""" ERROR: Could not find gdal-config, please ensure it is installed and on $PATH """)

    base_name = os.path.basename(dem_filename)
    base_name = os.path.splitext(base_name)[0]

    base_dir = user_output_dir + base_name + os.path.sep

    if not os.path.isdir(base_dir) and reuse_mesh:
        raise Exception(f"Using reuse_mesh=True, however the directory {base_dir} does not exist")

    # Delete previous dir (if exists)
    if os.path.isdir(base_dir) and not reuse_mesh:
        shutil.rmtree(base_dir, ignore_errors=True)



    # these have to be separate ifs for the logic to work correctly
    if not reuse_mesh:
        # make new output dir
        os.makedirs(base_dir)

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
        src_ds = gdal.Open(dem_filename)
        wkt = src_ds.GetProjection()
        srs = osr.SpatialReference()

        srs.ImportFromWkt(wkt)

        ext_str = ' -te %s %s %s %s -te_srs \"%s\" ' % (extent[0], extent[1], extent[2], extent[3], srs.ExportToProj4())
        src_ds = None

    cutline_cmd = ''
    if clip_to_shp:
        cutline_cmd = f' -cutline {clip_to_shp} -crop_to_cutline '

    e = 'GDAL_CACHEMAX=\"5%%\" %sgdalwarp %s %s -ot Float32 -multi -overwrite -dstnodata -9999 %s -t_srs \"%s\"' + ext_str
    subprocess.check_call([e % (gdal_prefix,
                                dem_filename, base_dir + base_name + output_file_name, cutline_cmd, srs_out.ExportToProj4())],
                          shell=True)

    if fill_holes:
        gt = src_ds.GetGeoTransform()

        pixel_width = gt[1]
        pixel_height = gt[5]

        dist = max([pixel_height, pixel_width]) * 5
        exec_str = '%sgdal_fillnodata.py %s -md %d' % (gdal_prefix, base_dir + base_name + output_file_name, dist)
        subprocess.check_call([exec_str], shell=True)

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

    exec_str = f'%sgdalwarp %s %s -ot Float32 -overwrite -multi -dstnodata -9999 -t_srs "%s" -te %s %s %s %s  -r '

    total_weights_param, use_weights_param = regularize_inputs(base_dir, exec_str, gdal_prefix, parameter_files,
                                                               pixel_height, pixel_width,
                                                               srs_out, xmax, xmin, ymax, ymin, fill_holes, nworkers_gdal)

    total_weights_ic, use_weights_ic = regularize_inputs(base_dir, exec_str, gdal_prefix, initial_conditions,
                                                         pixel_height, pixel_width,
                                                         srs_out, xmax, xmin, ymax, ymin, fill_holes, nworkers_gdal)

    use_weights = use_weights_param or use_weights_ic

    topo_weight = 1 - (total_weights_param + total_weights_ic)

    # over ride with user settings
    if user_no_weights:
        use_weights = False

    if use_weights and topo_weight < 0:
        raise RuntimeError("Parameter weights must equal 1")

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
            ingeom = feature.GetGeometryRef()
            geomBuffer = ingeom.Buffer(bufferDist)

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

    # if we aren't reusing the mesh, generate a new one
    if not reuse_mesh:
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

    # some paramters we want to use to constrain the mesh but don't actually want in the output. This let's use
    # remove them
    keys_to_drop = []
    for key, data in parameter_files.items():
        if 'drop' in data and data['drop'] == True:
            keys_to_drop.append(key)
    for key in keys_to_drop:
        print('Dropping parameter ' + key + ' and not writing to mesh due to user request')
        del parameter_files[key]

    # holds our main mesh structure which we will write out to json to read into CHM
    mesh = {'mesh': {}}
    mesh['mesh']['vertex'] = []

    read_header = False

    invalid_nodes = []  # any nodes that are outside of the domain.
    print('Reading nodes')
    with open(base_dir + 'PLGS' + base_name + '.1.node') as f:
        for line in f:
            if '#' not in line:
                if not read_header:
                    header = re.findall(r"[+-]?\d+(?:\.\d+)?", line)
                    num_nodes = int(header[0])
                    read_header = True
                    mesh['mesh']['nvertex'] = num_nodes
                else:
                    items = re.findall(r"[+-]?\d+(?:\.\d+)?", line)
                    mx = float(items[1])
                    my = float(items[2])

                    mz = extract_point(src_ds, mx, my)

                    if mz == dem.GetNoDataValue() or mz is None:
                        invalid_nodes.append(int(items[0]) - 1)

                    mesh['mesh']['vertex'].append([mx, my, mz])

    print('Length of invalid nodes = ' + str(len(invalid_nodes)))

    # read in the neighbour file, triangle topology
    print('Reading in neighbour file')
    read_header = False
    mesh['mesh']['neigh'] = []
    with open(base_dir + 'PLGS' + base_name + '.1.neigh') as elem:
        for line in elem:
            if '#' not in line:
                if not read_header:
                    read_header = True  # skip header, nothing to do here
                else:
                    items = re.findall(r"[+-]?\d+(?:\.\d+)?", line)
                    v0 = int(items[1]) - 1  # convert to zero indexing
                    v1 = int(items[2]) - 1
                    v2 = int(items[3]) - 1

                    mesh['mesh']['neigh'].append([v0, v1, v2])

    read_header = False

    print('Computing parameters and initial conditions')

    mesh['mesh']['elem'] = []
    mesh['mesh']['is_geographic'] = is_geographic

    # need to save the UTM coordinates so-as to be able to generate lat/long of points if needed later (e.g., CHM)
    if not is_geographic:
        mesh['mesh']['proj4'] = srs.ExportToProj4()
        mesh['mesh']['UTM_zone'] = srs.GetUTMZone()  # negative in southern hemisphere

    # holds parameters and initial conditions for CHM
    params = {}
    ics = {}
    for key, data in parameter_files.items():
        params[key] = []

    params['area'] = []
    params['id'] = []

    for key, data in initial_conditions.items():
        ics[key] = []

    # loop through all the triangles and assign the parameter and ic values to the triangle
    start_time = time.perf_counter()
    with open(base_dir + 'PLGS' + base_name + '.1.ele') as elem:
        for line in elem:
            if '#' not in line:
                if not read_header:
                    header = re.findall(r"[+-]?\d+(?:\.\d+)?", line)
                    nelem = int(header[0])
                    mesh['mesh']['nelem'] = nelem
                    read_header = True  # skip header
                else:
                    items = re.findall(r"[+-]?\d+(?:\.\d+)?", line)
                    v0 = int(items[1]) - 1  # convert to zero indexing
                    v1 = int(items[2]) - 1
                    v2 = int(items[3]) - 1

                    # if the node we have is invalid (out side of domain) we can try to fix it by interpolating from
                    # the surrounding nodes. this can happen when the outter domain constraint is effectively the
                    # limit of the DEM and the node ends up *just* outside of the domain due to numerical
                    # imprecision. estimate an invalid node's z coord from this triangles other nodes' z value
                    if v0 in invalid_nodes:
                        z_v1 = mesh['mesh']['vertex'][v1][2]
                        z_v2 = mesh['mesh']['vertex'][v2][2]
                        tmp = [x for x in [z_v1, z_v2] if x != dem.GetNoDataValue()]
                        # print 'found v0'
                        if len(tmp) != 0:
                            mesh['mesh']['vertex'][v0][2] = float(np.mean(tmp))
                            if verbose:
                                print('replaced invalid with ' + str(mesh['mesh']['vertex'][v0]))
                            invalid_nodes = [x for x in invalid_nodes if x != v0]  # remove from out invalid nodes list.

                    if v1 in invalid_nodes:
                        z_v0 = mesh['mesh']['vertex'][v0][2]
                        z_v2 = mesh['mesh']['vertex'][v2][2]
                        tmp = [x for x in [z_v0, z_v2] if x != dem.GetNoDataValue()]
                        # print 'found v1'
                        if len(tmp) != 0:
                            mesh['mesh']['vertex'][v1][2] = float(np.mean(tmp))
                            if verbose:
                                print('replaced invalid with ' + str(mesh['mesh']['vertex'][v1]))
                            invalid_nodes = [x for x in invalid_nodes if x != v1]  # remove from out invalid nodes list.

                    if v2 in invalid_nodes:
                        # print 'found v2'
                        z_v1 = mesh['mesh']['vertex'][v1][2]
                        z_v0 = mesh['mesh']['vertex'][v0][2]
                        tmp = [x for x in [z_v1, z_v0] if x != dem.GetNoDataValue()]
                        if len(tmp) != 0:
                            mesh['mesh']['vertex'][v2][2] = float(np.mean(tmp))
                            if verbose:
                                print('replaced invalid with ' + str(mesh['mesh']['vertex'][v2]))
                            invalid_nodes = [x for x in invalid_nodes if x != v2]  # remove from out invalid nodes list.
                    mesh['mesh']['elem'].append([v0, v1, v2])

    if len(invalid_nodes) > 0:
        errstr = 'Length of invalid nodes after correction= ' + str(len(invalid_nodes))
        errstr += 'This will have occurred if an entire triangle is outside of the domain. There is no way to ' \
                  'reconstruct this triangle. '
        errstr += 'Try reducing simplify_tol.'
        raise RuntimeError(errstr)

    gt = src_ds.GetGeoTransform()

    tris = range(mesh['mesh']['nelem'])

    ret_tri = []
    ret_tri_ic = []

    csize = len(tris) // nworkers
    if csize < 1:
        csize = 1

    with futures.ProcessPoolExecutor(max_workers=nworkers) as executor:
        for p, i in executor.map(
                partial(do_parameterize, gt, is_geographic, mesh, parameter_files, initial_conditions,
                        src_ds.RasterXSize, src_ds.RasterYSize, srs_out.ExportToProj4()), tris,
                chunksize=csize):
            ret_tri.append(p)
            ret_tri_ic.append(i)

    for key, data in parameter_files.items():
        parameter_files[key]['file'] = None

    for key, data in initial_conditions.items():
        initial_conditions[key]['file'] = None

    for t in ret_tri:
        for key, data in t.items():
            params[key].append(data)

    for t in ret_tri_ic:
        for key, data in t.items():
            ics[key].append(data)

    print('Total time took %s s' % str(round(time.perf_counter() - start_time, 2)))

    if output_write_vtu:
        print('Writing vtu file...')
        write_vtu(base_dir + base_name + '.vtu', mesh, params, ics)

    if output_write_shp:
        print('Writing shp file...')
        write_shp(base_dir + base_name + '_USM.shp', mesh, params, ics)

    print('Saving mesh to file ' + base_name + '.mesh')
    with open(user_output_dir + base_name + '.mesh', 'w') as outfile:
        json.dump(mesh, outfile, indent=4)

    print('Saving parameters to file ' + base_name + '.param')
    with open(user_output_dir + base_name + '.param', 'w') as outfile:
        json.dump(params, outfile, indent=4)

    print('Saving initial conditions  to file ' + base_name + '.ic')
    with open(user_output_dir + base_name + '.ic', 'w') as outfile:
        json.dump(ics, outfile, indent=4)
    print('Done')


def do_parameterize(gt, is_geographic, mesh, parameter_files, initial_conditions, RasterXSize, RasterYSize, srs_proj4,
                    elem):
    params = {}
    ics = {}
    srs_out = osr.SpatialReference()

    srs_out.ImportFromProj4(srs_proj4)

    for key, data in parameter_files.items():
        if data['file'] is None:
            parameter_files[key]['file'] = []
            for f in data['filename']:
                ds = gdal.Open(f)
                if ds is None:
                    raise RuntimeError('Error: Unable to open raster for: %s' % key)
                parameter_files[key]['file'].append(ds)

    for key, data in initial_conditions.items():
        if data['file'] is None:
            initial_conditions[key]['file'] = []
            for f in data['filename']:
                ds = gdal.Open(f)
                if ds is None:
                    raise RuntimeError('Error: Unable to open raster for: %s' % key)
                initial_conditions[key]['file'].append(ds)

    params['id'] = elem

    v0 = mesh['mesh']['elem'][elem][0]
    v1 = mesh['mesh']['elem'][elem][1]
    v2 = mesh['mesh']['elem'][elem][2]

    # Create a temporary vector layer in memory
    mem_drv = ogr.GetDriverByName('Memory')
    mem_ds = mem_drv.CreateDataSource('out')
    mem_layer = mem_ds.CreateLayer('poly', srs_out, ogr.wkbPolygon)

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

    feature = ogr.Feature(mem_layer.GetLayerDefn())
    feature.SetGeometry(tpoly)

    mem_layer.CreateFeature(feature)

    area = 0
    # if the output is geographic, we need to project to get a reasonable area
    if is_geographic:
        #use an equal area mollweide projection
        srs_moll = osr.SpatialReference()

        srs_moll.ImportFromProj4("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ")

        # go from what we are outputting (which is geographic) to moll to compute the area
        transform = osr.CoordinateTransformation(srs_out, srs_moll)
        p = tpoly.Clone()
        p.Transform(transform)

        area = p.GetArea()
    else:
        area = tpoly.GetArea()

    params['area'] = area

    # calculate new geotransform of the feature subset
    geom = feature.geometry()
    src_offset = bbox_to_pixel_offsets(gt, geom.GetEnvelope(), RasterXSize, RasterYSize)
    new_gt = (
        (gt[0] + (src_offset[0] * gt[1])),
        gt[1],
        0.0,
        (gt[3] + (src_offset[1] * gt[5])),
        0.0,
        gt[5]
    )

    # get the value under each triangle from each parameter file
    for key, data in parameter_files.items():
        output = []

        for f, m in zip(data['file'], data['method']):
            output.append(rasterize_elem(f, mem_layer, m, srs_out, new_gt, src_offset))

        if 'classifier' in data:
            fn = cloudpickle.loads(data['classifier'])
            output = fn(*output)
        else:
            output = output[0]  # flatten the list for the append below

        # we want to write actual NaN to vtu for better displaying
        if output == -9999:
            output = float('nan')

        if output is None and 'classifier' in data:
            print(f'Error: The user-supplied classifier function for {key} returned None which is not valid.')
            exit(1)

        params[key] = output

    for key, data in initial_conditions.items():
        output = []

        for f, m in zip(data['file'], data['method']):
            output.append(rasterize_elem(f, mem_layer, m, srs_out, new_gt, src_offset))

        if 'classifier' in data:
            output = cloudpickle.loads(data['classifier'])(*output)
        else:
            output = output[0]  # flatten the list for the append below

        if output == -9999:
            output = float('nan')

        if output is None and 'classifier' in data:
            print(f'Error: The user-supplied classifier function for {key} returned None which is not valid.')
            exit(1)

        ics[key] = output

    return params, ics


def regularize_inputs(base_dir, exec_str, gdal_prefix, input_files, pixel_height, pixel_width, srs_out,
                      xmax, xmin, ymax, ymin, fill_holes, max_workers):
    # ensure all the weights sum to 1
    total_weights = 0
    use_weights = False
    param_args = []
    ret = []
    for key, data in input_files.items():
        if 'weight' in data:
            total_weights += data['weight']
            use_weights = True
        if 'classifier' in data:
            input_files[key]['classifier'] = cloudpickle.dumps(input_files[key]['classifier'])

        param_args.append((base_dir, data, exec_str, gdal_prefix, key, pixel_height,
                           pixel_width, srs_out.ExportToProj4(), xmax, xmin, ymax, ymin, fill_holes))

    with futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for (r) in executor.map(_future_regularize_inputs, param_args):
            ret.append(r)
    for r in ret:
        key = r['key']
        input_files[key]['filename'] = r['filename']
        input_files[key]['file'] = None

        if not isinstance(input_files[key]['method'], list):
            input_files[key]['method'] = [input_files[key]['method']]

    return total_weights, use_weights


def _future_regularize_inputs(args):
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


# Print iterations progress
# http://stackoverflow.com/a/34325723/410074
def printProgress(iteration, total, prefix='', suffix='', decimals=2, barLength=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iterations  - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
    """
    filledLength = int(round(barLength * iteration / float(total)))
    percents = round(100.00 * (iteration / float(total)), decimals)
    bar = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")


# Zonal stats from here https://gist.github.com/perrygeo/5667173
def bbox_to_pixel_offsets(gt, bbox, rasterXsize, rasterYsize):
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1

    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1

    xsize = x2 - x1
    ysize = y2 - y1

    # only apply this correction if we are touching the underlying raster.
    if x1 < rasterXsize and y1 < rasterYsize:
        # deal with small out of bounds
        if x1 < 0:
            x1 = 0

        if y1 < 0:
            y1 = 0

        if x1 + xsize > rasterXsize:
            xsize = rasterXsize - x1

        if y1 + ysize > rasterYsize:
            ysize = rasterYsize - y1

    return x1, y1, xsize, ysize


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
    mz = rb.ReadAsArray(px, py, 1, 1)
    if mz is None:
        return rb.GetNoDataValue()
    mz = float(mz.flatten()[0])

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


def rasterize_elem(rds, mem_layer, aggMethod, srs, new_gt, src_offset):
    raster = rds.GetRasterBand(1)
    src_array = raster.ReadAsArray(*src_offset)

    # Rasterize it
    driver = gdal.GetDriverByName('MEM')
    mask = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
    mask.SetGeoTransform(new_gt)
    mask.SetProjection(srs.ExportToWkt())
    err = gdal.RasterizeLayer(mask, [1], mem_layer, burn_values=[1], options=['ALL_TOUCHED=TRUE'])
    if err != 0:
        raise Exception("Error rasterizing layer: %s" % err)

    # holds a mask of where the triangle is on the raster
    mask_arr = mask.ReadAsArray()

    # Mask the source data array with our current feature
    src_array[(mask_arr == 0) | (mask_arr == raster.GetNoDataValue())] = np.nan

    output = -9999.0

    if callable(aggMethod):
        output = float(aggMethod(src_array))
        if np.isnan(output):
            output = raster.GetNoDataValue()
    else:
        if aggMethod == 'mode':
            vals, counts = np.unique(src_array, return_counts=True)
            output = float(vals[np.argmax(counts)])

        elif aggMethod == 'mean':
            output = float(np.nanmean(src_array))
        elif aggMethod == 'max':
            output = float(np.nanmax(src_array))
        elif aggMethod == 'min':
            output = float(np.nanmin(src_array))
        else:
            print('\n\nError: unknown data aggregation method %s\n\n' % aggMethod)
            exit(-1)

    # feature.SetField(key, output)

    # testing code
    # if output < 0:
    #     print "Found < 0"
    #
    #     # testing code
    #     tri_id = str(feature.GetField('triangle'))
    #     print "Tri ID = " + tri_id
    #
    #     print masked
    #
    #     mem_drv = ogr.GetDriverByName('ESRI Shapefile')
    #     mem_ds = mem_drv.CreateDataSource(tri_id+'.shp')
    #     mem_layer = mem_ds.CreateLayer('poly', srs, ogr.wkbPolygon)
    #     mem_layer.CreateFeature(feature.Clone())
    #
    #     driver = gdal.GetDriverByName('GTiff')
    #     rvds = driver.Create(tri_id+'.tiff', src_offset[2], src_offset[3], 1, gdal.GDT_Float32)
    #     rvds.SetGeoTransform(new_gt)
    #     rvds.SetProjection(wkt)
    #     outband = rvds.GetRasterBand(1)
    #     outband.WriteArray(src_array)
    #     # gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1], options=['ALL_TOUCHED=TRUE'])
    #
    #     rvds = None
    #     mem_layer = None
    #     exit(1)

    return output


# To remove the dependency on gdal python bindings and thus the scripts, this is a slimmed down version of the
# gdal_polygonize.py code from
# https://github.com/OSGeo/gdal/blob/release/2.4/gdal/swig/python/scripts/gdal_polygonize.py
def gdal_polygonize(src_filename, mask, dst_filename):
    options = []
    src_band_n = 1
    src_ds = gdal.Open(src_filename)

    srcband = src_ds.GetRasterBand(src_band_n)

    mask_ds = gdal.Open(mask)
    maskband = mask_ds.GetRasterBand(1)

    dst_ds = ogr.Open(dst_filename, update=1)
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
        for key, data in parameter_files.items():
            output = data[elem]
            feature.SetField(key[0:10], output)

        for key, data in initial_conditions.items():
            output = data[elem]
            feature.SetField(key[0:10], output)

        layer.CreateFeature(feature)
    output_usm.FlushCache()
    output_usm = None  # close file


def write_vtu(fname, mesh, parameter_files, initial_conditions):
    vtu = vtk.vtkUnstructuredGrid()

    output_vtk = fname
    vtuwriter = vtk.vtkXMLUnstructuredGridWriter()
    vtuwriter.SetFileName(output_vtk)

    # check what version of vtk we are using so we can avoid the api conflict
    # http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Replacement_of_SetInput#Replacement_of_SetInput.28.29_with_SetInputData.28.29_and_SetInputConnection.28.29
    if vtk.vtkVersion.GetVTKMajorVersion() > 5:
        vtuwriter.SetInputData(vtu)
    else:
        vtuwriter.SetInput(vtu)

    vtu_points = vtk.vtkPoints()
    vtu_triangles = vtk.vtkCellArray()

    vtu_points.SetNumberOfPoints(len(mesh['mesh']['vertex']))

    vtu_cells = {'elevation': vtk.vtkFloatArray(),
                 'cellid': vtk.vtkFloatArray(),
                 'area': vtk.vtkFloatArray()
                 }
    vtu_cells['elevation'].SetName('elevation')
    vtu_cells['area'].SetName('area')
    vtu_cells['cellid'].SetName('cellid')

    for key, data in parameter_files.items():
        k = '[param] ' + key
        vtu_cells[k] = vtk.vtkFloatArray()
        vtu_cells[k].SetName(k)

    for key, data in initial_conditions.items():
        k = '[ic] ' + key
        vtu_cells[k] = vtk.vtkFloatArray()
        vtu_cells[k].SetName(k)

    for elem in range(mesh['mesh']['nelem']):
        v0 = mesh['mesh']['elem'][elem][0]
        v1 = mesh['mesh']['elem'][elem][1]
        v2 = mesh['mesh']['elem'][elem][2]

        vtu_points.SetPoint(v0, mesh['mesh']['vertex'][v0][0],
                            mesh['mesh']['vertex'][v0][1],
                            mesh['mesh']['vertex'][v0][2])

        vtu_points.SetPoint(v1, mesh['mesh']['vertex'][v1][0],
                            mesh['mesh']['vertex'][v1][1],
                            mesh['mesh']['vertex'][v1][2])

        vtu_points.SetPoint(v2, mesh['mesh']['vertex'][v2][0],
                            mesh['mesh']['vertex'][v2][1],
                            mesh['mesh']['vertex'][v2][2])

        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, v0)
        triangle.GetPointIds().SetId(1, v1)
        triangle.GetPointIds().SetId(2, v2)

        vtu_triangles.InsertNextCell(triangle)
        vtu_cells['elevation'].InsertNextTuple1((mesh['mesh']['vertex'][v0][2] +
                                                 mesh['mesh']['vertex'][v1][2] +
                                                 mesh['mesh']['vertex'][v2][2]) / 3.)

        vtu_cells['cellid'].InsertNextTuple1(elem)

        area = 1
        vtu_cells['area'].InsertNextTuple1(area)

        for key, data in parameter_files.items():
            output = data[elem]
            vtu_cells['[param] ' + key].InsertNextTuple1(output)

        for key, data in initial_conditions.items():
            output = data[elem]
            vtu_cells['[ic] ' + key].InsertNextTuple1(output)

    vtu.SetPoints(vtu_points)
    vtu.SetCells(vtk.VTK_TRIANGLE, vtu_triangles)
    for p in vtu_cells.values():
        vtu.GetCellData().AddArray(p)

    vtk_proj4 = vtk.vtkStringArray()
    vtk_proj4.SetNumberOfComponents(1)
    vtk_proj4.SetName("proj4")
    vtk_proj4.InsertNextValue(mesh['mesh']['proj4'])

    vtu.GetFieldData().AddArray(vtk_proj4)

    vtuwriter.Write()


if __name__ == "__main__":
    main()
