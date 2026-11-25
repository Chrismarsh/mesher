from osgeo import gdal


# Fill holes identified by the first band's GDAL mask in the staged raster. The search
# distance is measured in pixels on this grid.
def fill_nodata(path, distance):
    with gdal.Open(path, gdal.GA_Update) as ds:
        band = ds.GetRasterBand(1)
        gdal.FillNodata(band, band.GetMaskBand(), distance, 0)
