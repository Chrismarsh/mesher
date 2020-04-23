dem_filename = '../data/ideal_gaussianHill.tif'

max_area= 50000**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 5    # 5 m maxe RMSE between triangle and underlying elevation set to -1 to skip tolerance checks
min_area = 30**2     #triangle area below which we will no longer refine, regardless of max_tolerance

