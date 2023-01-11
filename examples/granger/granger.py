dem_filename = '../data/granger1m.tif'

max_area= 99999999**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 1    #1 m RMSE
min_area = 5**2     #triangle area below which we will no longer refine, regardless of max_tolerance

lloyd_itr=1
simplify=True
MPI_nworkers = 1