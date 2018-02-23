#!/bin/bash

function run
{
	cd $1
	python ../../mesher.py $1.py
	cd ..
}

run 'flat'
run 'flat_veg'
run 'gaussian_hill'
run 'ideal_ridge'
run 'ideal_ridge_low_tol'
run 'lloyd'
run 'uniform'
run 'flat_stream'
run 'flat_stream_dem'
run 'granger_low_veg_weight'
run 'granger_high_veg_weight'
run 'flow_accumulation'
run 'dem_smoothing'
run 'CAGEO_figure'