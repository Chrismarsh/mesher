
dem_filename = '../data/ideal_flat.tif'

max_area= 9999999999999999**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = -1    #1 -1 to skip tolerance checks
min_area = 100**2     #triangle area below which we will no longer refine, regardless of max_tolerance

lloyd_itr=100
MPI_nworkers = 1