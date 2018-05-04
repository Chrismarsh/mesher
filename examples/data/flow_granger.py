import richdem as rd
import numpy as np
import subprocess
from osgeo import gdal, ogr


dem = rd.LoadGDAL("dem.tif")
dem =  dem.astype(np.float32,copy=False)
#Fill depressions with epsilon gradient to ensure drainage
rd.FillDepressions(dem, epsilon = True, in_place=True)

#Get flow accumulation with no explicit weighting. The default will be 1.
accum = rd.FlowAccumulation(dem, method='Dinf')


rd.SaveGDAL('flow_accumulation_granger.tif',accum)

