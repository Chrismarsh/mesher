dem_filename='../data/chro_extent_lowRes.tif'

max_area=5000**2
max_tolerance=50
min_area=200**2

use_input_prj=False
lloyd_itr=100
simplify=True
simplify_tol=100
simplify_buffer=-50

parameter_files = {    
                     'flow_accumulation':{
                         'file':'../data/flow_accumulation.tif',
                         'method':'mean',
                         'tolerance':50
                     }
                   }


