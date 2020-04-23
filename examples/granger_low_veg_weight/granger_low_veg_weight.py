dem_filename = '../data/granger1m.tif'

max_area= 50000**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 5    #5m maxe RMSE between triangle and underlying elevation set to -1 to skip tolerance checks
min_area = 30**2     #triangle area below which we will no longer refine, regardless of max_tolerance

use_weights = True
topo_weight=0.8
weight_threshold = 0.8
parameter_files = {
    'landcover': {'file': '../data/eosd.tif',
                  'method': 'mode',
                  'weight':0.2,
                  'tolerance':.9}
}


lloyd_itr=1
simplify=True
