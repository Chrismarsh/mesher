import os
import sys
try:
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
except ModuleNotFoundError:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from mesher.mesher_utls.bootstrap_utils import ensure_mesher_on_path
ensure_mesher_on_path()
import cloudpickle
from mpi4py import MPI
import numpy as np
from osgeo import gdal, ogr, osr
import importlib.util
import importlib.machinery
from functools import lru_cache
from mesher.mesher_utls.mpi_progress import debug_log
from exactextract import exact_extract
from exactextract.feature import JSONFeatureSource
from exactextract.raster import GDALRasterSource

gdal.UseExceptions()  # Enable exception support
ogr.UseExceptions()
osr.UseExceptions()

def str2bool(s: str) -> bool:
    if s.lower() == 'true':
        return True

    return False

def bbox_to_pixel_offsets(gt, bbox, rasterXsize, rasterYsize):
    """Clip the geometry envelope to this raster's own affine pixel grid."""
    inverse = gdal.InvGeoTransform(gt)
    if inverse is None:
        raise ValueError('Cannot invert raster geotransform')
    pixels = [gdal.ApplyGeoTransform(inverse, x, y)
              for x in (bbox[0], bbox[1]) for y in (bbox[2], bbox[3])]
    xs, ys = np.asarray(pixels).T
    # Include the high-edge pixel as in the existing ALL_TOUCHED sampler.
    # Clamp both endpoints, not the origin of an already-sized window.
    x0 = max(0, min(rasterXsize, int(np.floor(xs.min()))))
    y0 = max(0, min(rasterYsize, int(np.floor(ys.min()))))
    x1 = max(0, min(rasterXsize, int(np.floor(xs.max())) + 1))
    y1 = max(0, min(rasterYsize, int(np.floor(ys.max())) + 1))
    return x0, y0, x1-x0, y1-y0


def _read_valid_samples(raster, src_offset):
    """Read one window with the same missing-data rules for both samplers."""
    # Integer class rasters cannot hold NaN. Mask raw samples before applying
    # scale/offset so the NoData sentinel cannot become an ordinary value.
    src_array = np.asarray(raster.ReadAsArray(*src_offset), dtype=np.float64).copy()
    valid = np.isfinite(src_array)
    nodata = raster.GetNoDataValue()
    if nodata is not None:
        valid &= src_array != nodata
    if raster.GetMaskFlags() != gdal.GMF_ALL_VALID:
        valid &= raster.GetMaskBand().ReadAsArray(*src_offset) != 0
    src_array[~valid] = np.nan

    scale, offset = raster.GetScale(), raster.GetOffset()
    src_array = src_array * (1.0 if scale is None else scale) + (0.0 if offset is None else offset)
    src_array[~np.isfinite(src_array)] = np.nan
    return src_array


class _ValidRasterSource(GDALRasterSource):
    """Give exactextract finite, scaled samples with raw NoData already masked."""

    def read_window(self, x0, y0, nx, ny):
        return np.ma.masked_invalid(_read_valid_samples(self.band, (x0, y0, nx, ny)), copy=False)

    def nodata_value(self):
        # The window mask is authoritative, including for scaled rasters.
        return None


def _parameter_value(key, values, classifier=None):
    """Classify one aggregate per raster, then normalize the output scalar."""
    output = classifier(*values) if classifier is not None else values[0]
    if output is None:
        raise ValueError(f'Error: The user-supplied classifier function for {key} returned None which is not valid.')
    return float('nan') if output == -9999 else float(output)


def rasterize_elem(rds, mem_layer, aggMethod, srs, new_gt, src_offset):

    if src_offset[2] <= 0 or src_offset[3] <= 0:
        return float('nan')
    src_array = _read_valid_samples(rds.GetRasterBand(1), src_offset)

    driver = gdal.GetDriverByName('MEM')
    mask = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
    mask.SetGeoTransform(new_gt)
    mask.SetProjection(srs.ExportToWkt())
    err = gdal.RasterizeLayer(mask, [1], mem_layer, burn_values=[1], options=['ALL_TOUCHED=TRUE'])
    if err != 0:
        raise RuntimeError('Error rasterizing layer: %s' % err)
    src_array[mask.ReadAsArray() == 0] = np.nan
    if not np.isfinite(src_array).any():
        return float('nan')

    if callable(aggMethod):
        output = float(aggMethod(src_array))
    elif aggMethod == 'mode':
        vals, counts = np.unique(src_array[np.isfinite(src_array)], return_counts=True)
        output = float(vals[np.argmax(counts)])
    elif aggMethod in ('mean', 'average'):
        output = float(np.nanmean(src_array))
    elif aggMethod == 'max':
        output = float(np.nanmax(src_array))
    elif aggMethod == 'min':
        output = float(np.nanmin(src_array))
    else:
        raise ValueError('Unknown data aggregation method %s' % aggMethod)

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


@lru_cache(maxsize=4)
def _load_config(configfile):
    # Load once per worker, including before unpickling callbacks referenced
    # through the config module. Re-executing it for every triangle is unsafe
    # for stateful configuration and needlessly expensive.
    loader = importlib.machinery.SourceFileLoader('config', configfile)
    module = importlib.util.module_from_spec(importlib.util.spec_from_loader('config', loader))
    sys.modules['config'] = module
    loader.exec_module(module)
    return module


def do_parameterize(gt, is_geographic, mesh,
                    parameter_files, current_parameter,
                    initial_conditions, RasterXSize, RasterYSize, srs_proj4,
                    elem, configfile):

    _load_config(configfile)

    params = {}
    ics = {}
    srs_out = osr.SpatialReference()

    srs_out.ImportFromProj4(srs_proj4)

    params['id'] = int(elem)

    v0 = mesh['mesh']['elem'][elem][0]
    v1 = mesh['mesh']['elem'][elem][1]
    v2 = mesh['mesh']['elem'][elem][2]

    # Create a temporary vector layer in memory
    mem_drv = ogr.GetDriverByName('MEM') or ogr.GetDriverByName('Memory')
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

    # get the value under each triangle from each parameter file
    # for key, data in parameter_files.items():

    key = current_parameter
    data = parameter_files[key]

    output = []

    if data.get('file') is None:
        raise RuntimeError(f'Parameter {key} has no open rasters in data[\"file\"]')
    methods = data['method'] if isinstance(data['method'], list) else [data['method']]
    if len(data['file']) != len(methods):
        raise ValueError(f'{key}: Specify one aggregation method per input raster')
    for f, m in zip(data['file'], methods):
        raster_gt = f.GetGeoTransform()
        src_offset = bbox_to_pixel_offsets(raster_gt, geom.GetEnvelope(), f.RasterXSize, f.RasterYSize)
        x, y = gdal.ApplyGeoTransform(raster_gt, src_offset[0], src_offset[1])
        new_gt = (x, raster_gt[1], raster_gt[2], y, raster_gt[4], raster_gt[5])
        output.append(rasterize_elem(f, mem_layer, m, srs_out, new_gt, src_offset))

    classifier = cloudpickle.loads(data['classifier']) if 'classifier' in data else None
    params[key] = _parameter_value(key, output or [float('nan')], classifier)

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

    config = _load_config(configfile)
    log_options = {'mesher_debug': getattr(config, 'mesher_debug', False)}
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
        stats = {'mode': 'majority', 'mean': 'mean', 'average': 'mean', 'min': 'min', 'max': 'max'}
        if method not in stats:
            raise ValueError(f'Unknown data aggregation method {method}')
        return stats[method]

    # Fast path: compute per-triangle stats once via exactextract. Parameters
    # with custom reducers retain the legacy ALL_TOUCHED sampler.
    # A classifier receives one aggregate per input raster, not individual pixels.
    # exactextract is much faster for zonal stats;
    if use_exactextract:
        features = build_tri_features(mesh)

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
            methods = data['method'] if isinstance(data['method'], list) else [data['method']]
            if any(callable(method) for method in methods):
                slow_params[key] = data
            else:
                exact_params[key] = data

        for key, data in exact_params.items():
            if len(data) == 0:
                continue
            debug_log(log_options, f'Rank {MPI.COMM_WORLD.rank} {key}: exactextract aggregation; '
                      'classifier, if supplied, runs once per triangle after aggregation')
            files = data['filename'] if isinstance(data['filename'], list) else [data['filename']]
            methods = data['method'] if isinstance(data['method'], list) else [data['method']]
            if not files or len(files) != len(methods):
                raise ValueError(f'{key}: Specify one aggregation method per input raster')
            classifier = cloudpickle.loads(data['classifier']) if 'classifier' in data else None
            values = []
            for fpath, method in zip(files, methods):
                stats = stats_for_method(method)
                debug_log(log_options, str(fpath))
                source = _ValidRasterSource(fpath)
                result = exact_extract(source, features, stats)
                # majority/min/max may return None for empty intersections.
                values.append([float('nan') if row['properties'].get(stats) is None
                               else float(row['properties'][stats]) for row in result])
                del source
            for tri_idx in range(mesh['mesh']['nelem']):
                ret_tri[tri_idx][key] = _parameter_value(
                    key, [raster_values[tri_idx] for raster_values in values], classifier)

        # Custom reducers still receive the legacy unweighted pixel array.
        for key, data in slow_params.items():
            debug_log(log_options, f'Rank {MPI.COMM_WORLD.rank} {key}: legacy ALL_TOUCHED aggregation; '
                      'classifier, if supplied, runs once per triangle after aggregation')
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

        debug_log(log_options, f'Rank {MPI.COMM_WORLD.rank} writing output pickle')
        with open(f'pickled_param_args_rets_{MPI.COMM_WORLD.rank}.pickle', 'wb') as f:
            cloudpickle.dump(ret_tri, f)
        os.remove(pickle_file)
        if disconnect:
            comm = MPI.Comm.Get_parent()
            comm.Disconnect()
        return

    for key, data in parameter_files.items():

        debug_log(log_options, f'Rank {MPI.COMM_WORLD.rank} {key}')
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

    debug_log(log_options, f'Rank {MPI.COMM_WORLD.rank} writing output pickle')
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
