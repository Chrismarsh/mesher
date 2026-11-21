Configuration
--------------

Mesher depends heavily upon GDAL to handle the geospatial data and the GDAL python bindings. Mesher's input rasters can be in any 1-band raster that GDAL can open. Input rasters need to have a coordinate system.

.. note:: 
   Only the items in the `Required`_ section are needed for a valid configuration file. However, most of the power of mesher requires use of the other features.

Required
=========
Configuration files are stored in a python file and passed as an argument to ``mesher`` on the command line. For example:
::

   python mesher example_mesher_config.py


Therefore this configuration file must be compliant python code, but as such can contain arbitrary python code.
Throughout the configuration file, areas and tolerances are taken to be in metres (m).


A minimum configuration file would be:
::

   dem_filename = '../data/ideal_flat.tif'

   max_area= 9999999**2  #Effectively unlimited upper area -- allow tolerance check to refine it further
   max_tolerance = 5    #1m max RMSE between triangle and underlying elevation
   min_area = 5**2     #triangle area below which we will no longer refine, regardless of max_tolerance
   errormetric = 'rmse'
   mesher_path = '../mesher'


.. confval:: dem_filename

   :type: path

The base DEM to use. Single band raster. The extent of ``dem_filename`` is used to define the area that is meshed. Input parameters are clipped to this extent. However, parameters need not cover the entire extent of the simulation domain. Must have a coordinate system. The coordinate system is used to define the coordinate system of the output mesh. This can be overridden with the ``wkt_out`` option -- see `Coordinate Systems`_.

This DEM can have missing values. In such a case, the domain to be meshed is determined as the largest continuous area. The use of missing values can be used to produce complex and non-rectangular areas. However, the more complex the domain, the greater the number of triangles required along the edge. See ``simplify`` on how to avoid this.

.. confval:: max_area

   :type: double

This is a constraint on the maximum size (m^2) of a triangle. Any triangle that is produce that is above this threshold will be rejected and refined. If you wish to set no upper limit on this, specify some very large number, such ``max_area = 9999999999999999``.

.. confval:: max_tolerance

   :type: double

The maximum difference (vertical distance) between the triangle and the underlying raster as measured by the ``errormeteric``. This can optionally be specified as ``-1`` to skip the tolerance checks, which can be useful in producing uniformly sized triangles.

.. image:: images/mesher_tolerances.jpg

.. confval:: min_area

   :type: double
   :default: ``pixel_width * pixel_height``

The area (m^2) at or below which a well-shaped triangle is accepted without
sampling or checking elevation, land-cover or other raster criteria, including
weighted criteria. The minimum-angle requirement still applies. For example,
with ``min_area = 3``, a well-shaped 2 m^2 triangle is accepted even if its raster
error would exceed the requested tolerance; a poorly shaped triangle of the
same area still requires refinement.

A good starting value is the area of a DEM cell. This is a stopping threshold
for raster-driven refinement, not a lower bound on triangle size: shape
refinement and constraint splitting can produce smaller triangles.

.. note::
   Although there is a default value for this, it really should be specified with a value that makes sense.

.. confval:: errormetric

   :type: string

Assigned an string value that determines the error metric to use.

'mean_tol' = Mean elevation difference 
'max_tol' = Max elevation difference
'rmse'  = RMSE tolerance 

The RMSE produces the best distribution of triangle sizes and does not penalized individual raster cells like ``max_tol`` does. RMSE should generally be used.

.. confval:: mesher_path

   :type: string

Set the path to the backend mesher executable. Optionally use the MESHER_EXE environment variable. This is almost certainly not something that needs to be set by the user and is for development purposes.

.. confval:: nworkers

    :type: int
    :default: Number of physical CPUs

Mesher will use the number of physical CPUs to parallel various tasks. This allows for fine-tuning that parameter.

.. confval:: nworkers_gdal

    :type: int
    :default: Number of physical CPUs

The initial step where all input parameters and initial conditions are made resampled and clipped to the input DEM
extent can be memory heavy. Therefore this limits, separate from the number of workers for the main parameterize loop,
how many parallel workers should be used to limit memory consumption.

Environment variables
=====================
.. confval:: MESHER_EXE
   
   :type: string

Instead of specifying a ``mesher_path`` in the configuration file, the environment variable ``MESHER_EXE`` may be set to the binary. If both ``MESHER_EXE`` and ``mesher_path`` are defined, the configuration file path takes precedence.

.. confval:: MESHER_NWORKERS

    :type: int

Specifies number of CPUs to use for parallel tasks. This will override a configuration file ``nworkers`` value.

.. confval:: MESHER_NWORKERS_GDAL

    :type: int

Specifies number of CPUs to use for parallel gdal tasks. This will override a configuration file ``nworkers_gdal`` value.

.. confval:: GDAL_CACHEMAX

    :type: string

GDAL's raster I/O functionality will, by default, cache reads and writes by up-to 5% of total memory to increase efficiency.
When dealing with large rasters and parallel processing, this can result in out-of-memory errors. For example, 32 processes
can cache 160% of total available RAM. By default, mesher will internally set the parallel processing cache to be no more
than ~30% of total RAM shared across all processes. This is about ~1% per process. For most systems this gives a good
blend of performance and memory use. However, if out-of-memory situations occur, try setting this lower or to 0. More information
can be found `on the gdalwarp documentation <https://trac.osgeo.org/gdal/wiki/UserDocs/GdalWarp#WarpandCacheMemory:TechnicalDetails>`_
and `wiki page <https://trac.osgeo.org/gdal/wiki/ConfigOptions#GDAL_CACHEMAX>`_.

Lloyd iterations
================

.. confval:: lloyd_itr

   :type: int

Enables *n* Lloyd iterations for serial meshing. MPI cavity remeshing requires ``lloyd_itr=0``. As per the `CGAL documentation <https://doc.cgal.org/latest/Mesh_2/index.html#secMesh_2_optimization>`_: the goal of this mesh optimization is to improve the angles inside the mesh, and make them as close as possible to 60 degrees. 100 iterations is a suggested amount. However, please note this does invalidate the numerical guarantees about the minimum amount of error introduced to the mesh.


.. |no_lloyd| image:: images/no_lloyd.jpg
   :width: 100%
   :align: middle
.. |lloyd| image:: images/lloyd.jpg
   :width: 100%
   :align: middle

+-------------+
|  No Lloyd   | 
+-------------+
| |no_lloyd|  | 
+-------------+
|   Lloyd     | 
+-------------+
|   |lloyd|   | 
+-------------+

.. confval:: verbose
   
   :type: boolean
   :default: False

Enables verbose output

.. confval:: mesher_debug

   :type: boolean
   :default: False

Enables detailed MPI meshing diagnostics: backend heartbeats, rejected patch
attempts, local cavity growth, preserved-work counts, clipping extensions and
per-merge timings. Also enables the existing C++ backend debug option.

By default MPI tile generation and each merge stage show compact completion
counts and elapsed time. Terminals update one progress line; redirected/PBS
logs receive start/end lines and a heartbeat every 30 seconds. Counts represent
completed tasks, not triangle counts or estimated remaining time. Fatal errors
and reproduction-folder paths remain visible with debug disabled. Failed
backend attempts retain their inputs and full ``backend.log`` in either mode.

Outputs
========
.. confval:: user_output_dir

   :type: path

Sets a user-defined output directory instead of the automatically generated folder name.

.. confval:: reuse_mesh 

   :type: boolean
   :default: False

 If a mesh was already generated, and only applying a new parametrization is required, enabling this skips the mesh generation step.

.. confval:: reuse_projected_parameters

   :type: boolean
   :default: False

Set ``reuse_projected_parameters = True`` to reuse rasters from ``parameter_files`` that
Mesher has already aligned to the DEM grid. For example, leave ``reuse_mesh = False``
when changing mesh tolerances to generate a new mesh from the same rasters. The DEM
is still processed and triangle parameters are calculated for the current mesh.

Mesher writes the rasters to ``projected_parameters/`` inside the DEM output directory,
including on runs with reuse disabled. Missing rasters are built as needed. Reuse
requires matching source files, target grid, GDAL version, resampling algorithm and
hole-filling settings. Keep the original inputs accessible so Mesher can check the paths,
sizes and modification times of the files GDAL reports. Sources without local file
timestamps are processed on every run. Older runs with randomly named projected
TIFFs need one preprocessing pass.

.. confval:: write_shp
   
   :type: boolean
   :default: True

Writes a .shp file corresponding to the produced triangulation. Is an expensive operation and can slow down production of very large meshes with many parameters.

.. confval:: write_vtu

   :type: boolean
   :default: True

Write a .vtu file that can be view in Paraview. 


Coordinate systems
===================

By default mesher uses the coordinate system of the DEM as the output mesh's coordinate system. 

.. confval:: use_input_prj 

   :type: boolean
   :default: True

Use the input file's projection. This is useful for preserving a UTM input. 

.. confval:: wkt_out
   
   :type: string

The output coordinate system can be set using this variable. This needs to be in well known text (wkt) format, e.g.,
::

    wkt_out = "PROJCS[\"North_America_Albers_Equal_Area_Conic\"," \
              "     GEOGCS[\"GCS_North_American_1983\"," \
              "         DATUM[\"North_American_Datum_1983\"," \
              "             SPHEROID[\"GRS_1980\",6378137,298.257222101]]," \
              "         PRIMEM[\"Greenwich\",0]," \
              "         UNIT[\"Degree\",0.017453292519943295]]," \
              "     PROJECTION[\"Albers_Conic_Equal_Area\"]," \
              "     PARAMETER[\"False_Easting\",0]," \
              "     PARAMETER[\"False_Northing\",0]," \
              "     PARAMETER[\"longitude_of_center\",-96]," \
              "     PARAMETER[\"Standard_Parallel_1\",20]," \
              "     PARAMETER[\"Standard_Parallel_2\",60]," \
              "     PARAMETER[\"latitude_of_center\",40]," \
              "     UNIT[\"Meter\",1]," \
              "     AUTHORITY[\"EPSG\",\"102008\"]]"



Domain simplification
======================

.. confval:: simplify

   :type: boolean
   :default: False

As described in the `algorithm <algorithm.html>`_ section, the input DEM defines the area to be meshed. If no-data values are present, then the largest continuous area defines the area to be meshed. As a result, complex basin shapes will likely result in the creation of many triangles along the complex edges. This option can  be used to simplify the basin outline. Setting ``simplify`` to ``True`` requires setting a value for ``simplify_tol``.

Once the domain to be meshed is determined (and is represented by a polyline), this polyline is simplified so-as to have no more than ``simplify_tol`` meters of error. By default will enable ``simplify_buffer``. See ``no_simplify_buffer``.

Can only be used with a ``wkt_out`` that is a projected CRS.

.. confval:: simplify_tol

   :type: double
   :default: 10 m

The maximum error (m) the polygon simplification of ``simplify`` can introduce. Be careful as too high a tolerance will cause triangles to be crated that are outside of the bounds of the raster.

.. confval:: simplify_buffer

   :type: double
   :default: -10 m

Sets a negative buffer (i.e., contracts the meshing domain) to give ``simplify_tol`` more room to work. That is, usage of the simplify tolerance without this will likely put triangles outside of the valid data domain. Using this allows for the simplification to result in triangles that exist within the data domain. You can disable the use of this buffer by setting ``no_simplify_buffer=True``. ``simplify_buffer`` is enabled by default when ``simplify=True`` is given.

.. confval:: no_simplify_buffer

   :type: boolean
   :default: False

Disables ``simplify_buffer`` when ``simplify=True``

.. confval:: extent

   :type: list[4]
   :default: []

A large DEM may be subset to ``[xmin ymin xmax ymax]``. These are given in the coordinate system of the input DEM.

.. confval:: clip_to_shp

   :type: str
   :default: None

A large DEM may be clipped to the geometry specified by a shape file. This should be a path to the shape file e.g., `/my/path/geom.shp'. The shp must be georeferenced and will be reprojected by gdalwarp as required.


Domain Repair
==============

.. confval:: fill_holes

    :type: boolean
    :default: False

A 1-pass of the GDAL `fill no data <https://gdal.org/programs/gdal_fillnodata.html>`_ algorithm. Maximum search distance is defined as ``max([pixel_height, pixel_width]) * 5``.

Input DEM smoothing
====================

If the DEM quality is poor or if triangles close to the elevation raster cell size is required. If the ``min_area`` is approximately equal to the cell size of the raster and tolerance parameter ensures triangles of this size are being produced, then in complex terrain the stair stepping of the raster (due to non-continuous first derivative; i.e., slope) impacts the mesh quality as shown below.

.. image:: images/dem_no_smoothing.jpg

Repeated smoothing iterations can be done. Each smoothing iteration using cubic spline and resamples by ``iter * smoothing_scaling_factor``.

 This is the result of 1 smoothing iteration. 

.. image:: images/dem_smoothing_2.jpg

This is the result of 2 smoothing iterations. 

.. image:: images/dem_smoothing_2.jpg

.. confval:: do_smoothing 

   :type: boolean
   :default: False

Smooths the input DEM. 

.. confval:: smoothing_scaling_factor 

   :type: double
   :default: 2.0

 Smoothing factor for above option.

.. confval:: max_smooth_iter 

   :type: int
   :default: 1

 Number of iterations to smooth over. Each iteration the smoothing magnitude increases by ``iteration * smoothing_scaling_factor``.


Parameters
===========

Parameters, such as vegetation cover, flow accumulation, or soils, can be used in mesher in two ways: 1) have their values assigned to a triangle or 2) also constrain the mesh generation, such as shown with the vegetation heights shown on the `algorithm <algorithm.html>`_ page.

No constraint
*************

Parameters are given by key-value pairs in a dictionary, where ``value`` is a dict that contains a file name and an aggregation method:
::

   parameter_files = { 'param_name' : {'file':'file.tif','method':'mean'}}


For example:
::

   parameter_files ={ 
      'soils' : {'file':'/path/to/soils.tif','method':'mean'},
      'canopy_height' : {'file':'/path/to/veg.tif','method':'mean'}
      }


.. confval:: param_name

   :type: string

This is the name of the parameter as it will appear in the final mesh output. I.e., what you want this to be called for use in a model.

.. confval:: file
   
   :type: string

This is a fully qualified path to the raster file


.. confval:: method

   :type: string
   :default: 'mean'

The method controls how the various raster cells that make up the triangle are combined together to give. One of ``mode`` or ``mean`` (``average`` works too). The ``mean`` is the mean of all values, whereas ``mode`` takes the value that has the greatest number of cells. The ``mode`` is good for classified data, such as soil or vegetation type.

Optionally, ``method`` may be a user-specified function that accepts a numpy array and returns one value.


.. confval:: do_cell_resize
   
   :type: bool
   :default: true

Mesher will resize the cells of a raster to match the input DEM's resolution. However this can result in producing a lot of small triangles around a coarse raster's cells if a tolerance for this raster is set. For example, if using a coarse vegetation class with a high resoution DEM, it can produce a lot of edge triangles around the vegetation classes. Therefore, if a tolerance on a parameter layer is set, then ``do_cell_resize=False`` is automatically set.

.. confval:: drop

    :type: bool
    :default: false

Some parameters may be used to constrain the mesh during generation but a user may not want it written to file. Setting ``drop=True`` will not write this parameter to file.

Classifier
**********

A user-specified function can be given as a classified function. This allows, after the ``method`` function has combined all the rasters into a single value for the triangle, to re-value this value. This function takes 1 value, and returns 1 value. It is called **after** the ``method`` function is called.

Each input raster is aggregated independently, then the classifier is called
once per triangle with those aggregates. With ``use_exactextract=True`` (the
default), built-in methods use exactextract even when a classifier is supplied:
``mean``/``average`` weight pixels by their covered fraction and ``mode`` selects
the class with the greatest total coverage. ``min``/``max`` use intersecting
valid pixels. On the projected grids used for MPI meshing, coverage weighting
corresponds to weighting by intersected area. Classification still happens
after aggregation; it is not applied to individual pixels.

If any input uses a callable ``method``, all inputs for that parameter retain
the legacy ``ALL_TOUCHED`` sampler and its unweighted arrays. Setting
``use_exactextract=False`` selects that sampler for every parameter. Routine
sampling messages identify the selected backend only with ``mesher_debug=True``.

Raster NoData values, invalid mask pixels and nonfinite samples are excluded
before aggregation and before applying the band's scale/offset. Valid zero
values are retained. An empty raster intersection or an aggregation with no
valid samples returns ``NaN``. Classifiers receive this missing value and should
handle it explicitly instead of mapping it through a catch-all class. For
example, ``if not np.isfinite(value): return np.nan`` preserves missing data.
For multiple inputs, a classifier can decide which missing inputs matter (for
example, a valid water mask may be sufficient even when canopy height is missing).

For example, perhaps a vegetation density metric as derived from remote sensing is to be converted to a canopy height. This could be done as:
::

   def Tree_cover_2_VegHeight(value):
       if value >= 50:
           value = 10 # set the tree hight to 10 m
       else:
           value = 0.1 # otherwise set low density to open
       return value

   parameter_files ={'CanopyHeight': {'file':'60N_120W_treecover2010_v3.tif',
                                       'method':'mean',
                                       'classifier':Tree_cover_2_VegHeight
                                       }
                     }


Alternatively, a binary tree/no-tree parameter could be derived
::

   def Tree_cover_2_Simple_Canopy(value):
       if value >= 50:
           value = 0 # not open
       else:
           value = 1 # open
       return value

  parameter_files ={ 'landcover': {'file':'60N_120W_treecover2010_v3.tif',
                                    'method':'mean',
                                    'classifier':Tree_cover_2_Simple_Canopy
                                  } }



Multiple input rasters can be combined in to a single parameter using a more complex classifier. This works by passing ``file`` and ``method`` a list of length *n* of the files and aggregation methods. The classifier then takes *n* arguments and returns a singlge value. The following examples shows using a water mask and a landcover map to make a water/open/treed dataset.
::

   def make_landcover(water,veg):
     if water == 1:
       return 0 #keep it water

     if veg < 0.2:
       return 1 #this will be clearing/open

     return 2 # this is tree cover

   parameter_files = { 'landcover':{'file':['waterMask.tif','Simard_Pinto_3DGlobalVeg_L3C.tif'],
                                    'method':['mode','mean'],
                                    'classifier':make_landcover}}



Using as constraint
*******************

By setting an optional ``tolerance`` value in the dictionary, a raster can be used to constrain the mesh generation. If ``method`` is ``mode``, then this is a fractional percent of the dominate cell values to cover the triangle area. Otherwise, it is RMSE in the units of the raster's value.

For example:
::

   parameter_files ={ 
      'soils' : {'file':'/path/to/soils.tif','method':'mean','tolerance':0.6},
      'canopy_height' : {'file':'/path/to/veg.tif','method':'mean','tolerance':2}
      }

Assuming the soils raster is a classification map, then each triangle must have 60% of one soil type **as well as** be within 2 m RMSE to the canopy height.

Any number of rasters may have a tolerance. Further, if used with a user-defined classifier, then the tolerance check occurs after the classifier has run

.. note::
   As more tolerances are added, or tolerances become tighter, more and more triangles will be produced. Past a certain point, it does not become meaningful to use an approximating mesh!  


Initial conditions
===================

Initial conditions work exactly the same way as paramters. However, they are written to a seperate file than the parameters and are intended to be used as intial conditions in a numerical model. They may be used to constrain the mesh.

::

   initial_conditions = { ... }


Shape file constraints
=======================

Shape files may be used to further constrain the mesh, for example to rivers or basin outlines. The line segments in the shape file are treated as barries that triangles cannot cross; thus triangle edges represent the shapefile edges exactly. It may be benificial to simplify these edges somewhat so-as to avoid the creation of many small triangles.

.. image:: images/constraint_shpfile.jpg

This is further shown in :ref:`examples:flat_stream`.

::

   constraints = { 'river_network' :
                  {
                     'file': '../data/Stream.shp'
                     'simplify':1 # will be in shp file's original units
                  }
            }


.. confval:: constraints.keyname.simplify

   :type: double

Amount to simplify the shapefile edges by. Measured as maximum error between old and new lines. In the units of the shp file.


MPI
======

By default Mesher will use MPI to launch the Mesher backend tasks.

.. confval:: MPI_nworkers

    :type: int
    :default: Number of cores on machine (e.g., 4, 8, 10, etc)

    Number of MPI workers. Each worker processes one tile or one merge at a time.
    For an external launcher, set its rank count to the same value.

.. confval:: mpi_num_tiles

    :type: int or None
    :default: None

    Exact number of spatial tiles for a fresh MPI mesh. Mesher factors this count
    into rows and columns to make tiles as close to square as the projected
    domain permits. A prime count gives a single row or column. With the
    default, the existing grid selection uses the available worker count.

    For example, this meshes 96 tiles with at most eight active tile workers.
    Each worker finishes one tile before starting its next. Completed tiles are
    stored on disk and pairwise merge jobs use the same worker limit.

    .. code-block:: python

        mpi_mesh = True
        MPI_nworkers = 8
        mpi_num_tiles = 96

    Increase the tile count to reduce the initial area per meshing task.
    More tiles add seams, merge jobs and intermediate files. Peak job memory
    also includes DEM preprocessing, the active seam geometry and raster caches.

.. confval:: mpi_parameter_batch_size

    :type: int
    :default: 32768

    Maximum number of faces a parameter worker processes per batch. Set this
    independently of the spatial tile count, including when using
    ``reuse_mesh=True``. Smaller batches can reduce per-worker memory during
    raster aggregation at the cost of more batch setup.

If Mesher is used on a cluster to process a large domain, the use of a job scheduler, such as SLURM, may be optimal.

.. confval:: MPI_exec_str

    :type: string
    :default: None

    Set this to a command to use to invoked the MPI job. For example
    `MPI_exec_str='./submit_job.sh job.sh'`
    where `submit_job.sh` invokes the queue submission, e.g.,
    `sbatch "$@"`
    and `job.sh` contains
    `srun --label --unbuffered  python "$@"`

    Parameterization appends the worker script, a shared input-manifest path,
    ``False`` and the configuration path to this command, quoting these paths.
    The launcher must wait until all workers finish. Its rank count must match
    ``MPI_nworkers``.

    Workers use one MPI launch and read shared vertex/element NumPy arrays
    through read-only memory maps. Each rank owns a contiguous face range;
    there is no separate preparation launch or full vertex-array pickle per
    rank. Temporary inputs and rank results live in a unique directory under
    the output directory, which must be accessible to every compute node.

    If `MPI_exec_str` is provided `MPI_nworkers` must also be provided.

.. confval:: mpi_mesh

    :type: bool
    :default: False

    Enable MPI tile meshing. Produces a single stitched mesh.

.. confval:: mpi_merge_method

    :type: str
    :default: "remesh"

    The only supported value is ``remesh``. It generates overlapping tiles and replaces artificial boundaries
    with validated local CGAL meshes. Requires a rebuilt C++ executable,
    Shapely >= 2, a projected working CRS and ``lloyd_itr=0``.
    The shared-edge implementation has been retired.
    Remove ``mpi_merge_snap_tol``, ``mpi_shared_edge_spacing`` and
    ``mpi_shared_edge_constraints`` from old MPI configurations and regenerate
    shared-edge intermediates with ``reuse_mesh=False``. See :doc:`mpi` for
    recovery and validation rules.

.. confval:: mpi_remesh_band_width

    :type: double or None
    :default: None

    Initial halo and removal-band width in projected coordinate units. The
    automatic value is the smaller of ``2*sqrt(max_area)`` and 20% of the
    smallest tile dimension. The area/coordinate units must be consistent.
    Retained faces affected by interface refinement or adjacent quality failures
    are absorbed locally on retry. Consecutive local retries double the number
    of neighbouring rings removed (1, 2, 4, then capped at 4), within retained components.
    Other recoverable failures double the
    current removal band. The physical domain and quality thresholds are unchanged.

.. confval:: mpi_remesh_max_attempts

    :type: int
    :default: 12

    Maximum number of attempts within a loaded seam region, including final
    interface recovery. The same limit applies separately to the number of
    seam-region attempts, including the first solve. A value of 1 allows a
    single solve at each level. Exhaustion leaves any previous output intact.

.. confval:: mpi_remesh_max_vertices

    :type: int
    :default: 500000

    CGAL vertex budget per tile or cavity invocation. Refinement fails when
    the budget is reached. This does not cap memory used by retained tile
    arrays, graph construction, or GEOS validation. See the separate active-face
    budget below.

.. confval:: mpi_remesh_max_active_faces

    :type: int
    :default: 2000000

    Limit on faces loaded into an in-memory seam region. Checked before loading
    input faces and between local repair attempts. If expansion exceeds this
    limit, the merge fails without publishing; increase it only when memory
    permits. This is not a byte limit: native refinement, candidate arrays,
    boundary geometry and GEOS validation also allocate memory.

The distributed MPI smoother has also been retired because it did not preserve
validated mesh criteria. Remove ``mpi_lloyd_boundary_tol`` and any nonzero
``mpi_global_lloyd`` from MPI configurations. ``mpi_global_lloyd=0`` is accepted
as a compatibility no-op. Serial ``lloyd_itr`` remains available.

.. confval:: use_exactextract

    :type: bool
    :default: True

    Use the exactextract Python API for parameters with built-in aggregation
    methods, including parameters with a classifier. Supported methods are
    ``mean`` (or ``average``), ``mode`` (coverage-weighted majority), ``min`` and
    ``max``. Unknown methods raise an error instead of silently computing a mean.

    Classifiers receive one aggregate per raster in the configured input order,
    once per triangle. Missing aggregates are passed as ``NaN``. If any input
    uses a callable ``method``, all inputs for that parameter use the legacy
    ``ALL_TOUCHED`` sampler, preserving the reducer's unweighted pixel arrays.
    Set to False to use that sampler for all parameters. Its unweighted pixel
    statistics can differ from exactextract's coverage-weighted statistics,
    including mode tie handling. Every sampler uses each parameter raster's
    own grid and NoData metadata.

    Each rank extracts all requested built-in statistics for a given raster in
    one pass and reuses those aggregates across parameters. Numeric results are
    released after their last consumer; nothing is cached across runs.
    Projected triangle areas are calculated from XY coordinates in bounded
    NumPy batches. Geographic meshes retain the existing equal-area projection
    calculation; elevations do not affect the horizontal area parameter.

    Earlier versions always used the legacy sampler for classifier parameters.
    With the default setting, those parameters now use coverage weighting;
    values can change near class boundaries and for small triangles. Existing
    classifier functions and their argument order do not need to change.
    Set ``reuse_mesh=True`` to regenerate parameters on a completed mesh
    without remeshing. This sampling change does not alter mesh refinement
    criteria or require rebuilding the C++ backend.

.. confval:: write_param_per_file

    :type: bool
    :default: False

    Write each parameter to its own JSON file (`<output>.param.<name>.json`)
    instead of a single `<output>.param` file. This reduces memory pressure
    for very large meshes.

Timing output
-------------

Mesher reports preprocessing, meshing (or mesh reuse), mesh preparation,
parameterization and export separately. Mesh preparation includes topology
validation and elevation sampling. Parameterization also reports shared-input
writing, MPI worker execution (including launch overhead) and result collection.
Those three times are included in the overall parameterization duration.

``Total run took`` covers configuration loading through completed exports.
The older ``Total time took`` message measured only parameterization and should
not be compared directly with the new total-run duration. Compare it with
``Parameterization took`` instead.

For a reproducible local worker benchmark against the pre-optimization code,
run ``PYTHONPATH=src python tests/benchmark_parameterization.py --repeats 3``
from a source checkout. ``--case batched-stats`` exercises different statistics
from one raster; the default shares one mean across three classifiers. The
benchmark verifies identical outputs but does not measure MPI startup, input
preparation, full-domain meshing or production-filesystem performance.

Recorded medians over three runs on 20,000 synthetic triangles were 8.27 s to
1.78 s when three classifiers shared a raster mean (4.65x), and 11.42 s to
2.85 s when requesting mean/min/max together (4.01x). Outputs matched exactly
in both cases. Raw timings are in ``parameterization-benchmark.json``; these
single-worker measurements are not estimates of whole-domain MPI speedup.

MPI pitfalls and guards
-----------------------

With ``mpi_num_tiles=None``, tiling balances worker utilization with physical tile shape using the
projected domain width and height. It does not round the rank count down to a
perfect square. For a square domain, 95 ranks produce a 9-by-10 grid and 96
ranks produce an 8-by-12 grid. Some allocations can still leave ranks unused
when filling them would create long thin tiles. This is a geometric heuristic;
it does not predict terrain-dependent refinement cost.

Before launching MPI work, the smallest tile edge must exceed 1.25 times the
remesh band width. If an explicit band is too wide, reduce the tile count
or ``mpi_remesh_band_width``. Small positive-area domain fragments are retained;
they are not discarded to improve utilization.

A completed run writes ``<base_name>_mpi_layout.json`` with its chosen grid and
final mesh filename. Reuse reads that manifest, so the number of workers can
change when regenerating parameters. Older completed runs without a manifest
still use the original square-grid naming rule and require their original
``MPI_nworkers`` value.












Disk-backed preparation and export
----------------------------------

Mesh preparation stores coordinates, connectivity and global neighbours as NPY
memory maps in a temporary directory beneath the mesh output directory. Existing
NPZ meshes (including compressed archives) are extracted incrementally, so
``reuse_mesh=True`` does not require loading their arrays into the Python heap.
Serial Triangle files are read row by row into the same representation.

Parameter workers process up to ``mpi_parameter_batch_size`` faces per batch
(32,768 by default) and write disjoint
ranges of shared parameter columns. The parent validates rank completion and
IDs without gathering triangle dictionaries. JSON retains the CHM schema;
VTU uses one piece with streamed binary arrays and global point IDs. Shapefile
export remains a face-at-a-time writer. These changes are automatic.

The shared output filesystem must have room for mesh arrays, parameter columns,
and temporary topology indexes. Memory-mapped pages are reclaimable file cache;
RSS can include those pages. Native patch meshing, raster caches and boundary
geometry still have their own memory requirements. Temporary preparation files
are removed when export completes or raises an exception.

Seam-region merging
-------------------

Each hierarchical merge, including the final merge, maps its input arrays from
disk. A strip around the seam and a four-ring triangle collar are loaded for
native refinement and parent repair. All other faces remain immutable on disk.
The cavity is the exact target domain minus their coverage. An outer-interface
check requires identical edge endpoints before joining the local result to the
interiors; no approximate welding or relaxed quality criteria are used.

If a disk edge has been split, Mesher first tries to subdivide its owning
triangle using the exact vertices from the seam. Each child must satisfy the
native quality criteria, and the new boundary must match the remaining disk
edges. When this repair fails, Mesher loads neighbouring faces and reuses the
completed local mesh. On the last expansion it first reserves space for the
usual collars. It then loads the connected disk interiors touching each collar
if they fit the remaining active-face budget. Unrelated components stay on
disk. The final mesh goes through the usual topology checks.

When a rejected local solve leaves no mesh to reuse, the final seam attempt
opens the remaining connected disk regions at the artificial cuts where they
fit the active-face budget. A small tile pair can therefore be loaded in full.

CGAL can finish its refinement queue while a face near a constrained edge
still fails the minimum angle. Mesher audits every interior face. For residual
shape failures it bisects the longest available unconstrained edge of each
face and reruns constrained refinement and the audit. Recovery continues
while the failure count decreases or, at equal count, the worst angle improves.
The vertex budget applies throughout. The quality report identifies a
``constrained_corner`` when two constraint rays enclose an in-domain angle
below the required minimum. Such a sector requires a change to a retained
interface or to the input geometry to meet that minimum.

Stream verification follows a connected chain of CGAL-constrained output
edges between the exact input vertices of each clipped stream segment.
Intermediate vertices may deviate from the segment by coordinate roundoff.
Missing links and disconnected endpoints fail validation.

Global topology checks and final compaction still make sequential passes over
the mesh, using disk-backed tables and bounded edge-sorting runs. NPZ output is
still readable by existing tools and old completed NPZ meshes can be reused.
This reduces heap memory at the cost of temporary disk space and additional IO;
it does not eliminate full-dataset disk reads/writes. Boundary geometry scales
with boundary complexity, and the loaded region scales with seam complexity.

Merge metadata reports ``disk_interior_faces``, ``peak_loaded_faces`` (input
region size, excluding subsequently generated candidates), ``seam_attempts``
and ``seam_work_reused`` (faces carried into expansion attempts, not unique final
faces). Existing refinement counters describe the final local repair, with
prior successful local reports under ``seam_attempt_reports``. Each report counts
``unmatched_disk_edges`` before disk-side repair. A successful subdivision
records ``disk_parents_subdivided``. Expansions record ``disk_expansion_faces``
and whether they closed affected disk components (``closing_disk_components``).

``tests/benchmark_streaming.py`` measures export heap allocations separately
from process RSS. In a synthetic 100,000 versus 1,000,000-face export, traced
Python/NumPy peak allocations were about 2.97 MB in both cases. This excludes
native meshing and is not a whole-pipeline memory or production timing claim.
