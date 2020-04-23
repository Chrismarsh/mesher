
dem_filename = '../data/ideal_flat.tif'

max_area= 9999999999999999**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 5    # 5m max RMSE between triangle and underlying elevation set to -1 to skip tolerance checks
min_area = 5**2     #triangle area below which we will no longer refine, regardless of max_tolerance

constraints = { 'river_network' :
					{
						'file': '../data/Stream.shp',
						'simplify':1 # will be in original projection units
					}
 			}
