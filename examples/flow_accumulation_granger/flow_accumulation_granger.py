dem_filename='../data/granger1m.tif'

max_area=5000**2
max_tolerance=10
min_area=5**2

lloyd_itr=100
simplify=True
simplify_tol=100
simplify_buffer=-50

parameter_files = {    
                     'flow_accumulation':{
                         'file':'../data/flow_accumulation_granger.tif',
                         'method':'mean',
                         'tolerance':500
                     }
                   }


