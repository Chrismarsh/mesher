MPI
-----

Mesher can use MPI to regularize the input rasters in parallel,
and to do the regridding of the parameter rasters to the triangles. Mesher can also do a type of parallel mesh
generation. A few considerations about this are described below.


Mesher will invoke the MPI processes itself
++++++++++++++++++++++++++++++++++++++++++++++++

The main mesher application is single process and mesher will invoke the MPI sub-tasks itself. By default,
MPI Spawn will be used. This is only supported with OpenMPI and assumes that mesher is running in a context where
it makes sense to spawn MPI ranks (e.g., HPC job).  In constrained HPC context, it may be better to have the main
mesher process running in a small job and spawn MPI sub-jobs in dedicated jobs. This is where ``MPI_exec_str`` is used.
``MPI_exec_str`` can be a job submission script. Using the MPI spawn can sometimes not work, so setting
``MPI_exec_str=f'mpirun -n {MPI_nworkers} python '`` may be preferred.

Mesh seam
+++++++++++

There are two merge methods. ``mpi_merge_method="shared"`` is the compatibility
default described below. The opt-in ``"remesh"`` method removes artificial
seam constraints and is described in the next section.

The mesh generator is, by default, serial. For large meshes, this can be very slow. Enabling multiple MPI workers
will split the domain into subdomains and mesh each subdomain independently. Points are placed
along the seam, approximately spaced at the smallest allowed triangle size. This can be controlled with
``mpi_shared_edge_spacing``. Setting a lower density can ensure this boundary seam is not over reprsented, but the
triangles may not be strictly Delaunay while crossing the boundary. This can result in a regular pattern of small
triangles. This is more evident in areas with low topographic variability than in areas with high variability.

Independent refinement may add extra points to only one side of a shared edge.
After loading the merged MPI arrays (and after optional global Lloyd iterations),
Mesher reconciles these subdivisions before sampling elevations or face parameters.
It splits coarse triangles using existing boundary vertices, rebuilds global
neighbours, and requires each resulting seam segment to have a matching triangle
on the opposite side. It rejects unresolved hanging nodes, coincident boundary
vertices, degenerate/duplicate triangles, non-manifold edges, and gaps or
same-side overlaps along repaired edges. Coordinates are not moved; local area changes are bounded by the
XY detection tolerance (``1e-8`` in mesh coordinate units) and floating-point
roundoff. Reconciliation does not guarantee Delaunay/angle quality or check
arbitrary overlaps and coverage of the intended domain.

The log reports the number of split faces and added triangles. Triangle counts
and ordering can change, so parameter extraction and partitioning must use the
resulting mesh. Existing ``.param`` files and partition metadata must not be
attached to a mesh with repaired geometry.

Local cavity remeshing
++++++++++++++++++++++

Rebuild ``meshercxx`` after updating the source; this mode uses a new C++ domain
meshing entry point. Install the updated Python dependencies (including Shapely
>= 2). In the configuration, set:
::

   mpi_mesh = True
   mpi_merge_method = "remesh"
   reuse_mesh = False
   lloyd_itr = 0
   mpi_global_lloyd = 0

Use a projected working CRS with area and length units consistent with
``min_area``/``max_area``. Start in a separate output directory to retain older
results. A shared-boundary NPZ cannot be reused as a remesh result. Subsequent
``reuse_mesh=True`` runs may reuse a completed remesh result to regenerate
parameters; changing geometry criteria requires fresh tiles.

Each tile is meshed with a halo. At a merge, only whole triangles inside each
tile's nonoverlapping core and outside the removal band are retained. The
cavity is the intended polygonal domain **minus the coverage of those retained
triangles**. It is not a rectangle filled over existing triangles, and triangles
are never selected by centroid ownership. This set difference is central to
preventing overlapping strips at joins.

CGAL remeshes the cavity from its boundary and the original physical constraints,
using the same raster, categorical, weighting, area and shape criteria as normal
meshing. Old seam vertices inside the cavity are discarded. Domain marking
classifies connected regions separated by constrained edges. The region touching
infinity is always outside; bounded regions are classified using their largest
face. This avoids misclassifying an almost-collinear exterior triangle when its
centroid rounds onto the wrong side of a staircase boundary, while preserving
holes, disconnected components and regions enclosed by physical constraints.
The refinement queue uses a strict ordering: largest area violations first,
then the smallest minimum angle once the area bound is met. Intersecting
constraint lines are noded before insertion. This preserves the domain supplied
to the MPI stage; it does not repair upstream domain preprocessing choices.

Raster tolerance failures request mandatory CGAL refinement, like maximum-area
failures. An ordinary shape-refinement request can be dropped near protected
acute constraint clusters; using that status for raster error previously let
refinement finish while the final check still found raster violations. Weighted
failures are also preserved in CGAL's cached quality record. The minimum-area
stop and the existing numeric/categorical comparisons remain in force. These
criteria are shared with serial meshing, so both paths receive the correction.

The retained interface is a contract: every retained edge facing the cavity must
have an exact matching patch edge. If CGAL splits an interface, the candidate is
discarded and the affected retained faces plus neighbouring rings are absorbed
into the cavity. The first local retry uses one ring; consecutive local retries
double the depth (1, 2, 4, 8, ...), stopping at the remaining component boundary.
This avoids moving an interface mismatch one row at a time until attempts run
out. Growth cannot cross an already removed cavity into another retained
component. Unresolved quality failures touching a retained interface trigger
the same local growth. This removes retained protrusions that leave narrow gaps
against a physical boundary without widening the entire seam band. A coordinate
roundoff allowance is used only to locate affected faces; vertices are never
moved or welded and all acceptance checks still apply. Other recoverable
failures widen the removal band. Rows and
columns use the same kernel; the later merge remeshes across earlier seams at
four-tile intersections. No old internal tile line is kept as a constraint.

Physical perimeter edges are classified using a roundoff-sized band around the
domain boundary. The **entire edge** must lie in that band; touching the boundary
at an endpoint is insufficient. This handles refined points on sloping physical
edges at large projected coordinates without incorrectly demanding a neighbour
outside the domain. Internal edges still require exact vertex/edge matching;
no coordinates are changed. The same classification applies to hole boundaries.

Before acceptance, validation checks numeric degeneracy, indexed topology,
positive-area overlaps (triangle area versus union area), symmetric difference
against the exact cavity, physical constraints, interface pairing and hanging
boundary nodes. Area comparisons allow scale-dependent floating-point roundoff.
The C++ backend re-evaluates its criteria on every output face. A candidate with
unresolved quality failures is rejected, including at sharp physical constraints
where those criteria cannot be met. The retry and vertex budgets make failures
explicit; there is no silent fallback to the shared-edge algorithm.

Backend or validation failures retain a ``cavity-*`` directory whose path is
reported in the exception. It contains ``domain.json``, ``constraints.json``,
``invocation.json`` (the command, working directory and tile/merge identity),
and the complete ``backend.log``. Successful calls remove their temporary
directory. Failed-attempt directories can be deleted after diagnosis. The log
is streamed to disk rather than accumulated in the MPI worker's memory.
When final quality checks fail, ``mesh.json.quality.json`` records each rejected
triangle's coordinates, area, minimum angle, constrained edges and failure reason
(area, angle or raster criterion), with counts by reason. No candidate mesh is
written. ``Remesh candidate rejected`` reports a failed attempt; it is not a
terminal merge failure unless recovery is exhausted or a vertex budget is hit.
Merge logs identify the tile/group and report the local neighbour-ring depth and
number of retained faces removed. The final rejected attempt does not announce
or perform growth that would never be retried. Attempt and vertex limits remain
hard limits; accelerated local growth does not guarantee every case converges.

The vertex-budget error reports the limit, current/input vertex counts and
refinement steps. It is distinct from a GDAL allocation or read error. Inspect
the log and boundary detail before increasing ``mpi_remesh_max_vertices``.
Triangle raster sampling clips windows against all four raster edges and
rejects failed reads before using their data. Disjoint windows supply no raster
samples; this does not extend DEM coverage or make out-of-domain geometry valid
for the final elevation export.

Only validated results are written, followed by metadata containing attempt,
retained-triangle and replacement-triangle counts. Parameters are extracted after
geometry is finalized. The current MPI smoothing implementation is disabled for
this path because it does not preserve the checked constraints and criteria.

This change addresses seam geometry and correctness. Pairwise stages still load
the two tile/group arrays, and the final export remains in memory. It does not
yet implement streaming assembly. Large cavities and complex boundaries can be
expensive; the logs expose retries rather than silently accepting a bad join.

Validation tests
~~~~~~~~~~~~~~~

Build the C++ binary, then run:
::

   PYTHONPATH=src:tests MESHER_TEST_EXE=/absolute/path/to/meshercxx \
       python -m unittest discover -s tests -v

The C++ tests cover overlapping tiles, a shifted seam on sloping terrain, rugged
terrain, holes, disconnected components, crossing physical constraints, a
four-tile junction, retry exhaustion and the vertex budget. With OpenMPI
available, a four-rank tile/merge/export test runs as well. Without
``MESHER_TEST_EXE``, C++ integration tests are explicitly skipped.

The C++ raster sampling and refinement-queue ordering regressions can also be built with
``cmake -S . -B build -DBUILD_TESTING=ON`` and ``cmake --build build``; run them
with ``ctest --test-dir build --output-on-failure``. They cover partial/disjoint
windows, rotated geotransforms, large coordinates, valid zeros, NoData,
failed source reads and the ordering requirements of CGAL's refinement queue.
The Python integration suite includes a reduced production staircase boundary
that previously refined exterior faces until the vertex budget was exhausted.
It also covers local recovery from angle failures near a retained protrusion,
preservation of unrelated retained faces, and rejection of the same geometry
when the protrusion is a physical domain boundary that cannot be removed.
Further regressions cover repeatedly subdivided interfaces, growth stopping at
existing cavities, and CGAL-refined sloping physical edges at large projected
coordinates. Genuine gaps and subdivided internal edges remain rejected.
A seven-vertex production-cavity reduction with a synthetic elevation anomaly
checks that raster refinement is not dropped at protected constraint clusters,
in both weighted and unweighted modes. C++ tests check mandatory numeric and
categorical refinement, cached weighted results, disabled numeric rasters and
the minimum-area stop.

Recovering an earlier MPI run
++++++++++++++++++++++++++++

If the previous run's intermediate directory and final stitched ``.npz`` are
still present, set the following in the original configuration and run Mesher:
::

   mpi_mesh = True
   reuse_mesh = True

Keep the original DEM/output paths and ``MPI_nworkers`` so the saved stitched
file can be located. Mesher skips tile triangulation, reconciles the loaded
arrays, and regenerates the exported mesh and face parameters. The stored
stitched NPZ is not modified. This workflow requires the original raster inputs
and configuration; it cannot recover a complete model from the ``.mesh`` alone.
Normal exported outputs are regenerated, so retain copies of any results you
need to keep. ``meshaudit --rebuild-neighbours`` remains a neighbour-only tool
and will refuse a mesh containing hanging nodes.

Zoom in of a seam

|image0-mpi|


|image1-mpi|

lloyd iterations
++++++++++++++++++
Global lloyd iterations breaks the alignment across the seams.

MPI mesh gen
+++++++++++++

To do the triangulation in parallel, ``mpi_mesh=True`` must be set.



.. |image0-mpi| image:: images/mpi_seam.png
.. |image1-mpi| image:: images/mpi_seam_large.png
