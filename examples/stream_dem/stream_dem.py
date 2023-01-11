

dem_filename = '../data/granger1m.tif'

max_area= 99999999**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 5    #5 m RMSE
min_area = 25**2     #triangle area below which we will no longer refine, regardless of max_tolerance


constraints = { 'river_network' :
					{
						'file': '../data/Stream.shp'
						# 'simplify':5 # will be in original projection units
					}
 			}

lloyd_itr=1
MPI_nworkers = 1