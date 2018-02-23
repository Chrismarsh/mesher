
#Example configuration file using the sample data

dem_filename = '../data/ideal_ridge.tif'

max_area= 50000**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 5    #1m maxe RMSE between triangle and underlying elevation set to -1 to skip tolerance checks
min_area = 30**2     #triangle area below which we will no longer refine, regardless of max_tolerance

use_weights = False
topo_weight=0.1
weight_threshold = 0.9
# parameter_files = {
#     'landcover': {'file': '../DEMs/eosd.tif',
#                   'method': 'mode',
#                   'tolerance':.9}
# }


# lloyd_itr=100
#Simplify the outter boundary allowing at most simplify_tol difference between boundary and simplified boundary
simplify=True
simplify_tol=500
simplify_buffer=-100
#If only the parameter files have changed, this will reuse the triangulation to reduce run time
reuse_mesh = False

#error metric for triangle tolerance. rmse is default
errormetric = 'rmse'

#path to the mesher binary
mesher_path = '../mesher'
