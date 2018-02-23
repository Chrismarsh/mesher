These data are used by the validation scripts in generating the meshes. 

- `ESOD.tif` is a Landsat classified vegetation dataset for a subset of Granger Creek in the Wolf Creek Reserach Basin, Yukon.
- `granger1m.tif` is the 1 m, LiDAR derived elevation for the Granger Creek subset.
- `makedata.py` generates idealized gaussian hill, ridge line, flat, etc. casses for the same area as `granger1m.tif`. This should run prior to running the the `ideal_` cases.

Use the `runall.sh` script to run all these tests cases.