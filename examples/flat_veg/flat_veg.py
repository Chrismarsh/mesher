
dem_filename = '../data/ideal_flat.tif'

max_area= 50000**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 5    # 5 m max RMSE between triangle and underlying elevation set to -1 to skip tolerance checks
min_area = 30**2     # triangle area below which we will no longer refine, regardless of max_tolerance

parameter_files = {
    'landcover': {'file': '../data/eosd.tif', # vegetation landcover
                  'method': 'mode',
                  'tolerance':.9}
}
MPI_nworkers = 1