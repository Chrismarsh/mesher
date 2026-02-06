dem_filename='../data/chro_extent_lowRes.tif'

max_area=5000**2
max_tolerance=50
min_area=200**2

use_input_prj=False
lloyd_itr=0
simplify=True
simplify_tol=100
simplify_buffer=-150

parameter_files = {    
                     'flow_accumulation':{
                         'file':'../data/flow_accumulation.tif',
                         'method':'mean',
                         'tolerance':50
                     }
                   }


MPI_nworkers=64
mpi_mesh=True
MPI_exec_str='mpirun -n 64 python '
nworkers_gdal=1
mpi_global_lloyd = 1


mpi_shared_edge_constraints=True
mpi_shared_edge_spacing = 200
