#!/bin/bash

function run
{
	cd $1
	mesher.py $1.py
	cd ..
}

run 'dem_smoothing'
run 'flat'
run 'flat_stream'
run 'flat_veg'
run 'flow_accumulation'
run 'flow_accumulation_granger'
run 'gaussian_hill'
run 'granger'
run 'granger_high_veg_weight'
run 'granger_low_veg_weight'
run 'ideal_ridge'
run 'ideal_ridge_low_tol'
run 'lloyd'
run 'stream_dem'
run 'uniform'
run 'CAGEO_figure'

