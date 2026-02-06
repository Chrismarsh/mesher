import os
import sys
try:
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
ensure_mesher_on_path()
import cloudpickle
import warnings
from mpi4py import MPI
import numpy as np
from osgeo import gdal, ogr, osr
import importlib
import json

gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()
osr.UseExceptions()

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

    scale = raster.GetScale() or 1.0
    offset = raster.GetOffset() or 0.0

    src_array = src_array * scale + offset

    output = -9999.0

    if callable(aggMethod):
        output = float(aggMethod(src_array))
        if np.isnan(output):
            output = raster.GetNoDataValue()
    else:
        if aggMethod == 'mode':
            vals, counts = np.unique(src_array[~np.isnan(src_array)], return_counts=True)
            if len(vals) == 0:
                output = raster.GetNoDataValue()
            else:
                output = float(vals[np.argmax(counts)])

        elif aggMethod == 'mean':
            output = float(np.nanmean(src_array))
        elif aggMethod == 'max':
            output = float(np.nanmax(src_array))
        elif aggMethod == 'min':
            output = float(np.nanmin(src_array))
        else:
            raise Exception('\n\nError: unknown data aggregation method %s\n\n' % aggMethod)


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


def do_parameterize(gt, is_geographic, mesh,
                    parameter_files, current_parameter,
                    initial_conditions, RasterXSize, RasterYSize, srs_proj4,
                    elem, configfile):

    X = importlib.machinery.SourceFileLoader('config', configfile)
    X = X.load_module()

    params = {}
    ics = {}
    srs_out = osr.SpatialReference()

    srs_out.ImportFromProj4(srs_proj4)

    params['id'] = int(elem)

    v0 = mesh['mesh']['elem'][elem][0]
    v1 = mesh['mesh']['elem'][elem][1]
    v2 = mesh['mesh']['elem'][elem][2]

    # Create a temporary vector layer in memory
    mem_drv = ogr.GetDriverByName('MEM')
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

    # if we are computing area, we must handle this slightly differ
    if current_parameter == 'area':
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

        params[current_parameter] = area
        return params  # exit early

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
    # for key, data in parameter_files.items():

    key = current_parameter
    data = parameter_files[key]

    output = []

    if data.get('file') is None:
        raise RuntimeError(f'Parameter {key} has no open rasters in data[\"file\"]')
    for f, m in zip(data['file'], data['method']):
        output.append(rasterize_elem(f, mem_layer, m, srs_out, new_gt, src_offset))

    if 'classifier' in data:
        fn = cloudpickle.loads(data['classifier'])
        output = fn(*output)
    else:
        # flatten the list for the append below
        output = output[0] if len(output) else float('nan')

    # we want to write actual NaN to vtu for better displaying
    if output == -9999:
        output = float('nan')

    if output is None and 'classifier' in data:
        raise Exception(f'Error: The user-supplied classifier function for {key} returned None which is not valid.')

    params[key] = float(output)

    # for key, data in initial_conditions.items():
    #     output = []
    #
    #     for f, m in zip(data['filename'], data['method']):
    #         output.append(rasterize_elem(f, mem_layer, m, srs_out, new_gt, src_offset))
    #
    #     if 'classifier' in data:
    #         output = cloudpickle.loads(data['classifier'])(*output)
    #     else:
    #         output = output[0]  # flatten the list for the append below
    #
    #     if output == -9999:
    #         output = float('nan')
    #
    #     if output is None and 'classifier' in data:
    #         print(f'Error: The user-supplied classifier function for {key} returned None which is not valid.')
    #         exit(1)
    #
    #     ics[key] = output

    return params #, ics

def main(pickle_file: str,
         disconnect: bool,
         configfile: str,
         *extra_args):
    # if called from SLURM, etc, these cli are coming in as strings
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    prepare_only = '--prepare-only' in extra_args

    if prepare_only:
        with open(pickle_file, 'rb') as f:
            prep = cloudpickle.load(f)

        verts = np.load(prep['verts_path'], mmap_mode='r')
        elems = np.load(prep['elems_path'], mmap_mode='r')
        ntris = len(elems)
        my_tris = np.array_split(np.arange(ntris), MPI.COMM_WORLD.size)
        tri_idx = my_tris[MPI.COMM_WORLD.rank]

        subset_mesh = {'mesh': {}}
        subset_mesh['mesh']['vertex'] = verts
        subset_mesh['mesh']['nvertex'] = len(verts)
        subset_mesh['mesh']['elem'] = elems[tri_idx]
        subset_mesh['mesh']['nelem'] = len(tri_idx)

        param_args = [
            prep['gt'],
            prep['is_geographic'],
            subset_mesh,
            prep['parameter_files'],
            prep['initial_conditions'],
            prep['RasterXSize'],
            prep['RasterYSize'],
            prep['srs_proj4'],
            prep.get('use_exactextract', True),
        ]

        out_pickle = f'pickled_param_args_{MPI.COMM_WORLD.rank}.pickle'
        with open(out_pickle, 'wb') as f:
            cloudpickle.dump(param_args, f)

        if disconnect:
            comm = MPI.Comm.Get_parent()
            comm.Disconnect()
        return

    # load our correct mesh subset file
    pickle_file = pickle_file.replace('RANK', str(MPI.COMM_WORLD.rank))

    with open(pickle_file, 'rb') as f:
        param_args = cloudpickle.load(f)

    gt, is_geographic, mesh, parameter_files, initial_conditions, RasterXSize, RasterYSize, srs_proj4, use_exactextract = param_args

    exactextract = None
    # exactextract is much faster for zonal stats; fall back if it's unavailable.
    if use_exactextract:
        try:
            from exactextract import exact_extract
            from exactextract.feature import JSONFeatureSource

            exactextract = exact_extract
        except Exception:
            warnings.warn('exactextract not available; falling back to per-pixel parameterization')
            use_exactextract = False

    ret_tri = [{} for _ in range(mesh['mesh']['nelem'])]

    def build_tri_features(mesh):
        features = []
        verts = mesh['mesh']['vertex']
        elems = mesh['mesh']['elem']
        for idx, tri in enumerate(elems):
            v0 = verts[tri[0]]
            v1 = verts[tri[1]]
            v2 = verts[tri[2]]
            coords = [[v0[0], v0[1]], [v1[0], v1[1]], [v2[0], v2[1]], [v0[0], v0[1]]]
            features.append({
                'type': 'Feature',
                'properties': {'id': idx},
                'geometry': {'type': 'Polygon', 'coordinates': [coords]}
            })
        return features

    def stats_for_method(method):
        if method == 'mode':
            return 'majority'
        return 'mean'

    # Fast path: compute per-triangle stats once via exactextract. Parameters
    # with classifiers fall back to the per-pixel path so semantics match.
    if use_exactextract and exactextract is not None:
        features = build_tri_features(mesh)
        if JSONFeatureSource is not None:
            srs = osr.SpatialReference()
            srs.ImportFromProj4(srs_proj4)
            srs_wkt = srs.ExportToWkt()

            features = JSONFeatureSource(features, srs_wkt=srs_wkt)

        ret_tri = [{} for _ in range(mesh['mesh']['nelem'])]

        exact_params = {}
        slow_params = {}
        for key, data in parameter_files.items():
            if key == 'area':
                continue
            if 'classifier' in data:
                slow_params[key] = data
            else:
                exact_params[key] = data

        for key, data in exact_params.items():
            if len(data) == 0:
                continue
            print(f'Rank {MPI.COMM_WORLD.rank} {key}')
            print(data)
            files = data['filename'] if isinstance(data['filename'], list) else [data['filename']]
            methods = data['method'] if isinstance(data['method'], list) else [data['method']]
            values = []
            for fpath, method in zip(files, methods):
                stats = stats_for_method(method)
                print(fpath)
                result = exactextract(fpath, features, stats)
                values.append([row.get('properties').get(stats) for row in result])
            for tri_idx in range(mesh['mesh']['nelem']):
                ret_tri[tri_idx][key] = values[0][tri_idx]

        # Slow path for classifier-based parameters.
        for key, data in slow_params.items():
            warnings.warn(f'Rank {MPI.COMM_WORLD.rank} {key} uses per-pixel classifier; '
                          f'exactextract is skipped for this parameter.')
            print(f'Rank {MPI.COMM_WORLD.rank} {key} (per-pixel classifier)')
            parameter_files[key]['file'] = []
            if key != 'area':
                files = data['filename'] if isinstance(data['filename'], list) else [data['filename']]
                for f in files:
                    ds = gdal.Open(f)
                    if ds is None:
                        raise RuntimeError(f'Error: Unable to open raster for: {key}')
                    parameter_files[key]['file'].append(ds)

            for elem in range(0, mesh['mesh']['nelem']):
                ret = do_parameterize(gt, is_geographic, mesh, parameter_files, key,
                                      initial_conditions, RasterXSize, RasterYSize,
                                      srs_proj4, elem, configfile)
                for k, d in ret.items():
                    ret_tri[ret['id']][k] = d
            parameter_files[key]['file'] = []

        for tri_idx in range(mesh['mesh']['nelem']):
            ret_tri[tri_idx]['id'] = tri_idx
            if 'area' not in ret_tri[tri_idx]:
                ret = do_parameterize(gt, is_geographic, mesh, parameter_files, 'area',
                                      initial_conditions, RasterXSize, RasterYSize,
                                      srs_proj4, tri_idx, configfile)
                ret_tri[tri_idx]['area'] = ret['area']

        print(f'Rank {MPI.COMM_WORLD.rank} writing output pickle')
        with open(f'pickled_param_args_rets_{MPI.COMM_WORLD.rank}.pickle', 'wb') as f:
            cloudpickle.dump(ret_tri, f)
        os.remove(pickle_file)
        if disconnect:
            comm = MPI.Comm.Get_parent()
            comm.Disconnect()
        return

    for key, data in parameter_files.items():

        print(f'Rank {MPI.COMM_WORLD.rank} {key}')
        parameter_files[key]['file'] = []

        # load the data
        if key != 'area':
            for f in data['filename']:
                ds = gdal.Open(f)
                if ds is None:
                    raise RuntimeError(f'Error: Unable to open raster for: {key}')

                parameter_files[key]['file'].append(ds)

        for elem in range(0, mesh['mesh']['nelem']):
            ret = do_parameterize(gt, is_geographic, mesh, parameter_files, key,
                                           initial_conditions, RasterXSize, RasterYSize,
                                           srs_proj4, elem, configfile)
            # ret looks
            # [{'id': 10354, 'area': 3485.748251695186, 'landcover': 4.0},
            #  {'id': 10355, 'area': 3740.141789605841, 'landcover': 4.0},
            #  {'id': 10356, 'area': 2721.806314367801, 'landcover': 4.0}]
            for k, d in ret.items():
                ret_tri[ret['id']][k] = d

            # ret_tri.append(ret)

        parameter_files[key]['file'] = []

    print(f'Rank {MPI.COMM_WORLD.rank} writing output pickle')
    # there is no way to return the uuid mangled filename + param name  so save it to a pickle which we can get later
    with open(f'pickled_param_args_rets_{MPI.COMM_WORLD.rank}.pickle', 'wb') as f:
        cloudpickle.dump(ret_tri, f)

    os.remove(pickle_file)

    # have been run from the MPI.spawn, so disconnect from parent
    if disconnect:
        comm = MPI.Comm.Get_parent()
        comm.Disconnect()


if __name__ == '__main__':
    main(*sys.argv[1:])
