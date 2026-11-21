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
from collections import Counter
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


def _projected_triangle_areas(mesh, chunk_size=65536):
    """Calculate XY areas in bounded batches, subtracting a local origin."""
    vertices = np.asarray(mesh['mesh']['vertex'], dtype=float)
    triangles = np.asarray(mesh['mesh']['elem'], dtype=np.intp)
    areas = np.empty(len(triangles), dtype=float)
    for start in range(0, len(triangles), chunk_size):
        points = vertices[triangles[start:start+chunk_size], :2]
        u = points[:, 1] - points[:, 0]
        v = points[:, 2] - points[:, 0]
        areas[start:start+chunk_size] = .5 * np.abs(u[:, 0]*v[:, 1] - u[:, 1]*v[:, 0])
    return areas


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

def rank_face_range(nfaces, rank, nranks):
    """Match array_split's contiguous ownership without a global index array."""
    if nranks <= 0 or not 0 <= rank < nranks:
        raise ValueError('Invalid parameterization rank layout')
    count, extra = divmod(nfaces, nranks)
    start = rank*count + min(rank, extra)
    return start, start + count + (rank < extra)


def _load_rank_mesh(prep, rank, nranks):
    vertices = np.load(prep['verts_path'], mmap_mode='r')
    elements = np.load(prep['elems_path'], mmap_mode='r')
    start, stop = rank_face_range(len(elements), rank, nranks)
    # Slice views stay file-backed. Only coordinates referenced by this rank's
    # faces need to be paged in; no full vertex array is copied into a pickle.
    mesh = {'mesh': {'vertex': vertices, 'elem': elements[start:stop],
                     'nvertex': len(vertices), 'nelem': stop-start}}
    return mesh, start, stop


def _write_worker_results(prep, ret_tri, start, stop):
    rank = MPI.COMM_WORLD.rank
    for index, triangle in enumerate(ret_tri, start):
        triangle['id'] = index
    output = {'rank': rank, 'start': start, 'stop': stop, 'triangles': ret_tri}
    with open(os.path.join(prep['output_dir'], f'result_{rank}.pickle'), 'wb') as stream:
        cloudpickle.dump(output, stream)


def main(pickle_file: str, disconnect: bool, configfile: str):
    if isinstance(disconnect, str):
        disconnect = str2bool(disconnect)

    config = _load_config(configfile)
    log_options = {'mesher_debug': getattr(config, 'mesher_debug', False)}
    with open(pickle_file, 'rb') as stream:
        prep = cloudpickle.load(stream)
    if prep['nranks'] != MPI.COMM_WORLD.size:
        raise ValueError(f"Parameterization expected {prep['nranks']} MPI ranks, "
                         f"but the launcher started {MPI.COMM_WORLD.size}")
    mesh, start, stop = _load_rank_mesh(prep, MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size)
    if 'result_paths' not in prep:
        rows = _parameterize_batch(prep, mesh, configfile, log_options)
        _write_worker_results(prep, rows, start, stop)
    else:
        from mesher.mesh_store import batches, BATCH_SIZE
        columns = {key: np.load(path, mmap_mode='r+') for key, path in prep['result_paths'].items()}
        for lo, hi in batches(stop-start, prep.get('batch_size', BATCH_SIZE)):
            local = {'mesh': dict(mesh['mesh'], elem=mesh['mesh']['elem'][lo:hi], nelem=hi-lo)}
            rows = _parameterize_batch(prep, local, configfile, log_options)
            for key, array in columns.items():
                array[start+lo:start+hi] = (np.arange(start+lo, start+hi) if key == 'id' else
                                           [row[key] for row in rows])
        for array in columns.values():
            array.flush()
        # A marker is published only after every column in this rank's range is complete.
        output = {'rank': MPI.COMM_WORLD.rank, 'start': start, 'stop': stop, 'complete': True}
        with open(os.path.join(prep['output_dir'], f'result_{MPI.COMM_WORLD.rank}.pickle'), 'wb') as stream:
            cloudpickle.dump(output, stream)
    if disconnect:
        MPI.Comm.Get_parent().Disconnect()


def _parameterize_batch(prep, mesh, configfile, log_options):
    gt, is_geographic = prep['gt'], prep['is_geographic']
    parameter_files, initial_conditions = prep['parameter_files'], prep['initial_conditions']
    RasterXSize, RasterYSize = prep['RasterXSize'], prep['RasterYSize']
    srs_proj4 = prep['srs_proj4']
    use_exactextract = prep.get('use_exactextract', True)

    ret_tri = [{} for _ in range(mesh['mesh']['nelem'])]
    projected_areas = None if is_geographic else _projected_triangle_areas(mesh)

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

        # Batch all requested statistics for each raster. Keep only numeric
        # columns, scoped to this worker invocation, until their last consumer.
        planned = {}
        requested = {}
        remaining = Counter()
        for key, data in exact_params.items():
            files = data['filename'] if isinstance(data['filename'], list) else [data['filename']]
            methods = data['method'] if isinstance(data['method'], list) else [data['method']]
            if not files or len(files) != len(methods):
                raise ValueError(f'{key}: Specify one aggregation method per input raster')
            inputs = [(os.fspath(path), stats_for_method(method)) for path, method in zip(files, methods)]
            planned[key] = inputs
            remaining.update(inputs)
            for path, stat in inputs:
                requested.setdefault(path, {})[stat] = None

        cached = {}
        for key, inputs in planned.items():
            data = exact_params[key]
            debug_log(log_options, f'Rank {MPI.COMM_WORLD.rank} {key}: exactextract aggregation; '
                      'classifier, if supplied, runs once per triangle after aggregation')
            classifier = cloudpickle.loads(data['classifier']) if 'classifier' in data else None
            values = []
            for fpath, stat in inputs:
                if (fpath, stat) not in cached:
                    debug_log(log_options, f'{fpath}: {", ".join(requested[fpath])}')
                    source = _ValidRasterSource(fpath)
                    result = exact_extract(source, features, list(requested[fpath]))
                    for name in requested[fpath]:
                        # majority/min/max may return None for empty intersections.
                        cached[fpath, name] = np.fromiter(
                            (float('nan') if row['properties'].get(name) is None
                             else row['properties'][name] for row in result),
                            dtype=float, count=mesh['mesh']['nelem'])
                    del result, source
                values.append(cached[fpath, stat])
            for tri_idx in range(mesh['mesh']['nelem']):
                ret_tri[tri_idx][key] = _parameter_value(
                    key, [float(raster_values[tri_idx]) for raster_values in values], classifier)
            del values
            for pair in inputs:
                remaining[pair] -= 1
                if not remaining[pair]:
                    del cached[pair]

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
                if projected_areas is not None:
                    ret_tri[tri_idx]['area'] = float(projected_areas[tri_idx])
                else:
                    ret = do_parameterize(gt, is_geographic, mesh, parameter_files, 'area',
                                          initial_conditions, RasterXSize, RasterYSize,
                                          srs_proj4, tri_idx, configfile)
                    ret_tri[tri_idx]['area'] = ret['area']

        return ret_tri

    for key, data in parameter_files.items():

        debug_log(log_options, f'Rank {MPI.COMM_WORLD.rank} {key}')
        if key == 'area' and projected_areas is not None:
            for elem, area in enumerate(projected_areas):
                ret_tri[elem].update(id=elem, area=float(area))
            continue
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

    # Always publish area and ID, including an empty parameter configuration.
    for index, row in enumerate(ret_tri):
        row['id'] = index
        if 'area' not in row:
            row['area'] = (float(projected_areas[index]) if projected_areas is not None else
                           do_parameterize(gt, is_geographic, mesh, parameter_files, 'area',
                                           initial_conditions, RasterXSize, RasterYSize, srs_proj4,
                                           index, configfile)['area'])
    return ret_tri


if __name__ == '__main__':
    main(*sys.argv[1:])
