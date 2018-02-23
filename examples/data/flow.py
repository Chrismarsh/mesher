import richdem as rd
import numpy as np
import subprocess
from osgeo import gdal, ogr
dem = rd.LoadGDAL("chro_extent_lowRes.tif")
# dem = rd.rdarray(dem, no_data=-9999)
#Fill depressions with epsilon gradient to ensure drainage
rd.FillDepressions(dem,  in_place=True)

#Get flow accumulation with no explicit weighting. The default will be 1.
accum_d8 = rd.FlowAccumulation(dem, method='Dinf')
# accum_d8[ accum_d8 < 5005] = 0
# accum_d8[ accum_d8 > 0] = 1
# d8_fig = rd.rdShow(accum_d8, zxmin=450, zxmax=550, zymin=550, zymax=450, figsize=(8,5.5), axes=False, cmap='jet')


rd.SaveGDAL('flow_accumulation.tif',accum_d8)


# tmp_raster = 'd8.tif'
# base_dir=''
# plgs_shp='rivernetwork.shp'
# subprocess.check_call(['gdal_polygonize.py %s -b 1 -mask %s -f "ESRI Shapefile" %s' % (tmp_raster, tmp_raster,
#                                                                                        base_dir +
#                                                                                        plgs_shp)], shell=True)


# exec_string = 'ogr2ogr -overwrite %s %s  -nlt LINESTRING' % ('line_' + plgs_shp, base_dir + plgs_shp)

# # if simplify:
# simplify_tol=50
# exec_string = exec_string + ' -simplify ' + str(simplify_tol)

# subprocess.check_call(exec_string, shell=True)


# exit(1)
# src_ds = gdal.Open( "d8.tif" )
# if src_ds is None:
#     print 'Unable to open %s' % src_filename
#     sys.exit(1)

# try:
#     srcband = src_ds.GetRasterBand(1)
# except RuntimeError, e:
#     # for example, try GetRasterBand(10)
#     print 'Band ( %i ) not found' % band_num
#     print e
#     sys.exit(1)

# dst_layername = "rivernetwork"
# drv = ogr.GetDriverByName("ESRI Shapefile")
# dst_ds = drv.CreateDataSource( dst_layername + ".shp" )
# dst_layer = dst_ds.CreateLayer(dst_layername, srs = None )

# gdal.Polygonize( srcband, None, dst_layer, -1, [], callback=None )
# dst_ds=None