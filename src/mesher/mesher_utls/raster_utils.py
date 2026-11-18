from osgeo import gdal


def fill_nodata(path, distance):
    with gdal.Open(path, gdal.GA_Update) as ds:
        band = ds.GetRasterBand(1)
        gdal.FillNodata(band, band.GetMaskBand(), distance, 0)
