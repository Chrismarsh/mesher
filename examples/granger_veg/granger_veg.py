
def make_landcover(land,water):
    outtype = land

    if water != 1 or land == 0:
        outtype = 18

    return outtype

def water(w):
    out = 1
    if w != 1:
        out = 0

    return out

dem_filename = '../data/granger1m.tif'

max_area= 50**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
max_tolerance = 0.25    # 5m maxe RMSE between triangle and underlying elevation set to -1 to skip tolerance checks
min_area = 1*2     #triangle area below which we will no longer refine, regardless of max_tolerance

use_weights = True
reuse_mesh = False

MPI_nworkers=4
mpi_mesh=True
MPI_exec_str='mpirun -n 4 python '
nworkers_gdal=1
mpi_global_lloyd = 0


mpi_shared_edge_constraints=True
mpi_shared_edge_spacing = 500


parameter_files = {
    'landcover': {'file': '../data/eosd.tif',
                  'method': 'mode'},

    # more complex workflows
    #'watermask' : {'file':'water_mask/Hansen_GFC-2019-v1.7-datamask.vrt', 'method':'mode', 'classifier':water, 'tolerance':0.43},
    # 'landcover': {'file': ['land_cover/canada_2015_v2/CAN_NALCMS_2015_v2_land_cover_30m/CAN_NALCMS_2015_v2_land_cover_30m.tif',
    #                    'water_mask/Hansen_GFC-2019-v1.7-datamask.vrt'], 'method':['mode','mode'],'classifier':make_landcover}

}


# lloyd_itr=
simplify=True
