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

MPI mesh generation
+++++++++++++++++++

Set ``mpi_mesh=True`` to generate overlapping tiles in parallel. Cavity
remeshing is the default and only supported merge method. It replaces
artificial tile seams while preserving validated triangles away from repairs.
The former shared-edge algorithm, its boundary densification and snap merge,
and its export-time triangle splitter have been retired.

Remove ``mpi_shared_edge_spacing``, ``mpi_shared_edge_constraints`` and
``mpi_merge_snap_tol`` from old MPI configurations. An explicit
``mpi_merge_method="shared"`` fails with a migration message; omit it or set
``"remesh"``. Old shared-edge intermediates cannot be reused: generate fresh
tiles with ``reuse_mesh=False``. Keep older results in a separate output
directory if needed.

Export builds global neighbours and checks topology before sampling elevations
or parameters. It preserves triangle count, connectivity and XY coordinates;
it rejects hanging nodes and coincident boundary vertices instead of trying to
repair the geometry after quality validation.

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

Use a projected working CRS with area and length units consistent with
``min_area``/``max_area``. Start in a separate output directory to retain older
results. A shared-boundary NPZ cannot be reused as a remesh result. Subsequent
``reuse_mesh=True`` runs may reuse a completed remesh result to regenerate
parameters; changing geometry criteria requires fresh tiles.

Each tile is meshed with a halo. Only connected portions of that halo which
intersect the tile core with positive area are meshed. A core with no domain
area produces an empty tile, even when its halo catches a boundary fragment.
This avoids unnecessary angle failures where a simplified/buffered boundary
clips the halo into an acute triangle. The global domain remains unchanged;
the owning tiles and subsequent merges cover those fragments. Point or line
contact does not count as ownership, and no small-area cutoff is applied.
Owned fragments still have to satisfy all shape requirements.

An owned component can also have an acute corner created by the rectangular
cut. Refinement cannot remove an input angle below the shape bound. At such
corners, extend the clipping window locally with a small polygonal disk so its
new edge meets the physical boundary at a wider angle. Intersect that enlarged
window with the unchanged physical domain. This preserves all owned core area,
holes and physical corners; it does not trim slivers or relax angle checks.
The disk radius is bounded by half the halo/band width and a quarter of either
incident edge length. At most four geometry passes are made before the normal
strict backend checks. Logs report the number of extensions and added halo area.

Intermediate merge groups use the same local extension where their core box
would recreate an acute cut. Their rectangular ownership boxes stay unchanged.
These temporary halo triangles are subject to normal whole-triangle ownership
at the next merge. The final group covers the complete physical domain, so no
artificial outer cuts or extra halo area remain. This clipping correction is
Python-only; it does not need a further C++ rebuild after the preserved-work
backend has been installed.

At a merge, only whole triangles inside each tile's nonoverlapping core and
outside the removal band are retained. The
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
failures are also preserved in CGAL's cached quality record. At or below
``min_area``, a well-shaped triangle is accepted before any raster sampling;
numeric, categorical and weighted criteria cannot force further refinement.
Shape quality is still required at every area. Above ``min_area``, the existing
numeric/categorical comparisons remain in force. These criteria are shared with
serial meshing, so both paths receive the correction.

The retained interface is a contract: every internal edge must have an exact
matching edge on the other side. CGAL's constraint hierarchy tracks each retained
edge by ID and returns its ordered subdivisions. When an edge splits,
``meshercxx`` first triangulates the affected parent boundaries in one batch,
without inserting additional points. Children must preserve the parent footprint
and pass the same quality criteria. Multiple split sides and multiple inserted
points are handled together; coordinates are never moved or welded.

If subdivision fails quality, remove the affected parents and adjacent patch
triangles plus one neighbouring ring. Keep validated generated triangles and
unaffected original triangles outside this region, and mesh only the missing
coverage on the next attempt. Later repairs operate on this working mesh.
Unchanged triangles keep their quality result within this merge, under the same
raster inputs/settings; there is no persistent quality cache between runs.
A patch whose own final quality check failed is rejected in full.

Promotion is transactional: validate subdivision footprints and prescribed
boundaries, then check the retained remainder before changing the private working
state. Incomplete working states are never published. Fully validated backend
patches receive one full coverage/overlap check, followed by exact interface and
hanging-node checks at assembly. Unresolved quality failures touching retained
interfaces still grow the local cavity, with depth capped at four rings for the
production backend. Other recoverable failures widen the band. The final
component recovery remains a last resort. Rows and columns use the same kernel;
later merges can repair earlier seams at four-tile intersections.

Before invoking CGAL, the merge checks retained interface vertices against all
constraint segments and constraint vertices against retained interface segments.
A nonincident point closer than the coordinate-roundoff tolerance can form a
valid GEOS polygon but an unresolved slit for CGAL's inexact constructions.
Such a candidate raises ``PatchPrecisionError`` and selects the implicated
retained owners for the existing local growth path. Exact incidences and
resolved gaps are left alone. This tolerance only selects whole triangles for
removal; it never snaps vertices, fills physical holes or relaxes acceptance
criteria. The physical domain and constraints stay unchanged.

Optimized CMake builds retain CGAL structural assertions (``CGAL_DEBUG``).
This prevents refinement from continuing with invalid handles after a failed
invariant: the patch entry point catches the exception and writes diagnostics.
This native safeguard requires rebuilding ``meshercxx``. Actual backend signal
termination remains fatal and is not treated as a recoverable geometry failure.

Local growth can keep exposing another edge that refinement splits. If a retry
is still needed before the final allowed attempt, all retained components with
a retained/cavity interface are absorbed into that last patch. This removes
every artificial interface rather than spending the last attempt on another
finite ring expansion. Complete disconnected regions bounded only by the
physical perimeter are preserved. The log reports ``Final interface recovery``
and the number of removed and retained faces. This can remesh most or all of the
current merge group and increase time and memory use. The same vertex budget,
raster/shape criteria, coverage and topology checks still apply. A budget
failure stops immediately; it never triggers a larger recovery patch.

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
``invocation.json`` (the command, working directory, tile/merge identity and core bounds),
``interfaces.json`` and the complete ``backend.log``.
For a precision rejection, ``precision.json`` records interface IDs and the
coordinate tolerance, and ``backend.log`` explicitly says the backend was not
started. The same reproduction inputs are retained for inspection.
Subdivision failures preserve a ``repair-*`` directory with parent boundaries,
invocation and log. ``mesh.json.stats.json`` records backend counters on both
success and ordinary failure. Successful calls remove their temporary
directory. Failed-attempt directories can be deleted after diagnosis. The log
is streamed to disk rather than accumulated in the MPI worker's memory.
When final quality checks fail, ``mesh.json.quality.json`` records each rejected
triangle's coordinates, area, minimum angle, constrained edges and failure reason
(area, angle or raster criterion), with counts by reason. No candidate mesh is
written. ``Remesh candidate rejected`` reports a failed attempt; it is not a
terminal merge failure unless recovery is exhausted, a vertex budget is hit, or
the backend is terminated by a signal. Signal termination stops retries immediately.
Merge logs identify the tile/group and report the local neighbour-ring depth and
number of retained faces removed. The final rejected attempt does not announce
or perform growth that would never be retried. Attempt and vertex limits remain
hard limits. The final interface recovery resolves the moving artificial
boundary, but physical constraints, quality requirements or resource limits can
still prevent a valid result.

The vertex-budget error reports the limit, current/input vertex counts and
refinement steps. It is distinct from a GDAL allocation or read error. Inspect
the log and boundary detail before increasing ``mpi_remesh_max_vertices``.
Triangle raster sampling clips windows against all four raster edges and
rejects failed reads before using their data. Disjoint windows supply no raster
samples; this does not extend DEM coverage or make out-of-domain geometry valid
for the final elevation export.

Only validated results are written. Metadata records attempts, original triangles
retained, generated triangles reused, successfully subdivided parents, generated
patch triangle count, local repair count and ``interface_recovery``.
``boundary_halo_area`` records temporary area outside an intermediate group's
core; it is zero for the completed whole-domain merge. Per-attempt
records include cavity area, preserved work and rejection reason. Timings cover
preparation, backend execution, refinement, raster evaluation, validation,
local repair, stitching and I/O. Backend/refinement/raster timings are nested;
local repair time includes its backend and validation. They must not all be added
as independent phases. Backend progress/heartbeats are printed about every ten
seconds while waiting; counters measure work, not a convergence percentage.
Parameters are extracted after geometry is finalized. The distributed MPI
smoother has been removed because it did not preserve the checked constraints
and criteria. Remove ``mpi_lloyd_boundary_tol`` and any nonzero
``mpi_global_lloyd`` from MPI configurations; ``mpi_global_lloyd=0`` remains a
compatibility no-op. Serial Lloyd smoothing is still available via ``lloyd_itr``.

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
four-tile junction, retry exhaustion and the vertex budget. Buffered boundary
tests reproduce both unowned three-vertex fragments and owned five-vertex acute
cuts, including at large projected coordinates. A complete tile/group merge
checks that local halo extensions preserve the domain and disappear from the
final result; actual acute physical corners still fail quality. With OpenMPI
available, a four- and six-rank tile/merge/export tests run as well, including a rectangular
grid with an odd number of columns. Without
``MESHER_TEST_EXE``, C++ integration tests are explicitly skipped.

The supplied 463-vertex ``stage0_r3_c1`` cavity is also a regression fixture.
It reproduces the reported angle failures without raster inputs, verifies that
recovery selects all three retained tips near the group boundary, and checks a
controlled local replacement with complete coverage and no overlap. This tests
the rejection/recovery path; a rejected candidate alone does not mean the MPI
run failed. Replaying the entire production retry also requires its retained
triangle connectivity and raster inputs.

The exact 288-vertex ``stage0_r5_c2`` input from ``cavity-55lyskn1`` is committed
as ``tests/fixtures/retained_near_touch_cavity.geojson``. Without rasters it
reproduces a release-build segfault in CGAL 6.0.1's edge refinement: a retained
tip is approximately 1.7e-10 metres from the sloping outer boundary. Preflight
selects its two incident interfaces before starting C++. A controlled removal
of that tip exposes three ordinary angle failures; removing those retained
tips then produces a strictly checked patch with complete coverage and no
overlap with the remaining retained footprint. A separate complete merge
test uses the same near-touch triangle and confirms that local recovery keeps
an unrelated original triangle and does not widen the band. A direct native
replay bypasses preflight to check that optimized builds report failure without
a signal. These regressions do not require production rasters, but do not
replay the complete production MPI run or its original retained connectivity.

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
With a short attempt limit, a moving-interface regression requires the final
component recovery and preserves an unchanged disconnected island. A real CGAL
regression forces interface subdivision and verifies recovery while preserving
a physical hole and stream constraint. Budget failures and invalid replacement
coverage still fail without modifying the inputs.
A seven-vertex production-cavity reduction with a synthetic elevation anomaly
checks that raster refinement is not dropped at protected constraint clusters,
in both weighted and unweighted modes. C++ tests check mandatory numeric and
categorical refinement, cached weighted results, disabled numeric rasters and
the minimum-area stop. An unreadable raster regression verifies that triangles
at or below the stopping threshold perform no raster reads, in both weighting
modes and for numeric/categorical inputs. Above the threshold, read errors still
fail explicitly; small triangles with poor angles remain rejected.

Recovering an earlier MPI run
++++++++++++++++++++++++++++

If the previous run's intermediate directory and final stitched ``.npz`` are
still present, set the following in the original configuration and run Mesher:
::

   mpi_mesh = True
   reuse_mesh = True

Keep the original DEM/output paths. New runs record their grid and final mesh
filename in ``<base_name>_mpi_layout.json``; reuse reads that completed layout
even if ``MPI_nworkers`` has changed. For older runs without this manifest,
keep the original ``MPI_nworkers`` so the saved stitched file can be located. Only a completed cavity-remesh run with matching metadata
is eligible. Mesher skips tile triangulation, validates the loaded topology,
and regenerates the exported mesh and face parameters. The stored
stitched NPZ is not modified. This workflow requires the original raster inputs
and configuration; it cannot recover a complete model from the ``.mesh`` alone.
Normal exported outputs are regenerated, so retain copies of any results you
need to keep. ``meshaudit --rebuild-neighbours`` remains a neighbour-only tool
and will refuse a mesh containing hanging nodes.

Preservation regressions cover multiple split sides, large coordinates, the
minimum-area stop, raster failures above that stop, poor-angle children, unchanged
patch interiors across retries, rollback after a later failure, invalid child
coverage and backend signal termination.

A reproducible synthetic benchmark compares two tile workers followed by one
pairwise merge with a historical implementation and a serial whole-domain patch
reference (not the complete serial CLI preprocessing/export pipeline):
::

   git show 23b9658:src/mesher/mesher_utls/mpi_remesh.py > /tmp/baseline_mpi_remesh.py
   PYTHONPATH=src python tests/benchmark_mpi_preservation.py \
       --backend /path/to/new/meshercxx \
       --baseline-backend /path/to/old/meshercxx \
       --baseline-module /tmp/baseline_mpi_remesh.py \
       --size 1200 --output /tmp/preservation-benchmark.json

Results include wall/CPU time, sampled Linux process-tree PSS, triangle counts
and repeated cavity work. This is a local regression benchmark, not evidence of
cluster scaling or production TVC speedup. Real DEM/Hansen inputs and retained
tile arrays are still needed for that comparison.

One local run (CGAL 6.0.1, GDAL 3.10.3, identical 1200 x 600 synthetic domain,
maximum area 8, minimum area 0.25, elevation tolerance 1) produced:

.. list-table:: Synthetic preservation benchmark
   :header-rows: 1

   * - Method
     - Generation wall time
     - Merge wall time
     - Merge attempts
     - Generated patch triangles
   * - Previous implementation (23b9658)
     - 38.77 s
     - 35.50 s
     - 11
     - 357,077
   * - Preserved work
     - 6.46 s
     - 3.36 s
     - 4
     - 926
   * - Serial whole-domain patch reference
     - 5.67 s
     - N/A
     - N/A
     - N/A

The preserved-work merge retained 135,687 of the initial 135,735 retained
triangles and reused 967 generated triangles in its final result. It subdivided
38 parents and did not use broad recovery. Generation was about six times faster
than the previous MPI strategy, while the serial reference remained faster in
this case. The final meshes have approximately 137,000 triangles and all pass
independent coverage/topology checks. Raw CPU, sampled process-tree PSS and
per-attempt measurements are in ``mpi-preservation-benchmark.json``. Timing is
machine-dependent and is not a production performance guarantee.
