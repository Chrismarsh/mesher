# mesher agent notes

## What this codebase does
Mesher is a hydrology-focused, multi-objective unstructured mesh generator. It takes a DEM plus optional
feature constraints (streams, landcover, etc.), preprocesses inputs with GDAL, and generates a variable
resolution triangular mesh with quality controls (e.g., Lloyd optimization and tolerance checks). Output
formats include shapefiles and VTK/VTU, with helper tools for conversions and mesh stats.

## Architecture at a glance
- Python CLI and pipeline: `src/mesher/mesher_cli.py` reads a user config file, prepares raster/vector
  inputs with GDAL, and drives the mesh workflow.
- C++ core: `src-cxx/` builds the `mesher` executable that performs the actual triangulation and
  tolerance checks (uses GDAL, CGAL, Boost).
- Utilities and tools: `src/mesher/tools/` contains converters and analyzers
  (`mesh2shp`, `mesh2vtu`, `meshmerge`, `meshpermutation`, `meshstats`).

## Key entry points
- CLI entry points are declared in `pyproject.toml` under `[project.scripts]`.
- Main user flow is `mesher` (Python), which expects a config file path.
- The config reader (`read_config` in `src/mesher/mesher_cli.py`) documents required options such as
  `dem_filename` and `max_area`, plus optional constraints and smoothing behavior.

## Build and environment notes
- Build uses `scikit-build-core` + CMake to compile the C++ core.
- Runtime dependencies include GDAL, VTK, CGAL, Boost, METIS, MPI (see `pyproject.toml` and README).
- There are example/test scripts in the repo root (e.g., `test.py`, `test-sub.py`) and docs in `docs/`.

