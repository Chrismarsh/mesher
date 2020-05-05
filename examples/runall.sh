#!/bin/bash
set -x
set -e
function run
{
	cd $1
	mesher.py $1.py
	cd ..
}

run 'dem_smoothing'
run 'flat'
if [[ -z "${TRAVIS}" ]]; then
  run 'flat_stream'
fi
run 'flat_veg'

#don't run these on travis as they take too long
if [[ -z "${TRAVIS}" ]]; then
  run 'flow_accumulation'
  run 'flow_accumulation_granger'
fi
run 'gaussian_hill'
run 'granger'
if [[ -z "${TRAVIS}" ]]; then
  run 'granger_high_veg_weight'
fi

run 'granger_low_veg_weight'
run 'ideal_ridge'
run 'ideal_ridge_low_tol'
run 'lloyd'
run 'stream_dem'
run 'uniform'
run 'CAGEO_figure'

