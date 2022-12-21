import os
import sys
import cloudpickle
from mpi4py import MPI
import numpy as np
import uuid
from osgeo import gdal, ogr, osr
import subprocess

def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True

    return False

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
            vals, counts = np.unique( src_array[~np.isnan(src_array)], return_counts=True)
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

def do_parameterize(args, elem):
    gt, is_geographic, mesh, parameter_files, initial_conditions, RasterXSize, RasterYSize, srs_proj4 = args

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

    return params #, ics

def main(pickle_file: str,
         disconnect: bool):
    # if called from SLURM, etc, these cli are coming in as strings
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    my_tris = np.array_split([x for x in range(param_args[2]['mesh']['nelem'])], MPI.COMM_WORLD.size)
    my_tris = my_tris[MPI.COMM_WORLD.rank]

    ret_tri = []

    for elem in my_tris:
        ret_tri.append(do_parameterize(param_args, elem))

    # there is no way to return the uuid mangled filename + param name  so save it to a pickle which we can get later
    with open(f'pickled_param_args_rets_{MPI.COMM_WORLD.rank}.pickle', 'wb') as f:
        cloudpickle.dump(ret_tri, f)

    # have been run from the MPI.spawn, so disconnect from parent
    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
