#!/usr/bin/env python
import numpy as np
import os

# Write a text file containing information on the WN maps used as input to mesher

# Directory where are stored to WN maps
dir_wind = '/home/vvionnet/data3/Run_WindNinja/WN_Lajoie/param_lajoie'

# Direcrtory where are stored mesher command
dir_mesh = '/home/vvionnet/snow_models/mesher/lajoie_srtm'

# Name of the wind maps:
name_run_WN = 'ref-DEM-utm'

# Number of wind maps to be used
ncat = 12

############################
#No change after these lines
############################

# Angle between each wind maps
delta_wind = 360/ncat

list_wind = np.arange(0,360,delta_wind)

# List of variables to be inclused in mesher
list_var = ['U','V','spd_up']

dic_suffix = {'U':'_U','V':'_V','spd_up':''}

# File containing results
fic_res = dir_mesh+"/config_WN.txt"

if(os.path.isfile(fic_res)):
    os.remove(fic_res)

file1 = open(fic_res,"w")

for ii,ww in enumerate(list_wind):

    if(ii==0):
       ii=ncat

    for var in list_var:
        ll = "'Ninja"+str(ii)+dic_suffix[var]+"' : {'file':'"+dir_wind+"/"+name_run_WN+"_"+str(ww)+"_"+var+".vrt','method':'mean'}, \n"
        file1.write(ll)
        
    file1.write("\n")
     
    
    
 

