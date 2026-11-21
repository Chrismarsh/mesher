"""Merge an expandable seam region; immutable interiors remain file-backed.

Full-mesh passes are sequential/batched IO (selection, topology, publication).
Only the seam region and its interface/physical-boundary geometry are assembled
in memory. The existing native refinement and parent-repair rules still apply.
"""
from contextlib import ExitStack
import math
import os
from pathlib import Path
import tempfile
import time
import zipfile

import numpy as np
from shapely import union_all
from shapely.geometry import LineString, box
from shapely.ops import polygonize

from mesher.mesh_store import batches, disk_array, open_mesh_arrays, check_boundary_geometry
from mesher.mesh_topology import build_neighbours
from . import mpi_remesh as remesh
from .mpi_progress import debug_log


def selected_boundary(vertices, triangles, neighbours, active):
    """Yield boundaries of the inactive part, including its active-side cuts."""
    for lo, hi in batches(len(triangles)):
        block = neighbours[lo:hi]
        retained = ~active[lo:hi]
        for slot in range(3):
            adjacent = block[:, slot]
            cut = adjacent < 0
            valid = ~cut
            cut[valid] = active[adjacent[valid]]
            for local in np.flatnonzero(retained & cut):
                face = lo + int(local)
                tri = triangles[face]
                p, q, r = vertices[tri, :2]
                i, j = tri[(slot+1) % 3], tri[(slot+2) % 3]
                a, b = vertices[i, :2], vertices[j, :2]
                if (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0]) < 0:
                    i, j, a, b = j, i, b, a
                yield face, (i, j), (tuple(a), tuple(b)), adjacent[local] >= 0


def inactive_coverage(vertices, triangles, neighbours, active):
    directed, edges, interfaces = set(), [], {}
    cut_vertices = {}
    for face, ids, pair, internal in selected_boundary(vertices, triangles, neighbours, active):
        directed.add(pair)
        edges.append(LineString(pair))
        if internal:
            interfaces[tuple(sorted(pair))] = face
            cut_vertices.update(zip(ids, pair))
    pieces = []
    for polygon in polygonize(edges):
        coords = list(polygon.exterior.coords)
        if not polygon.exterior.is_ccw:
            coords.reverse()
        if (coords[0], coords[1]) in directed:
            pieces.append(polygon)
    coverage = remesh.polygon_domain(union_all(pieces))
    total = 0.
    for lo, hi in batches(len(triangles)):
        pts = vertices[triangles[lo:hi][~active[lo:hi]], :2]
        u, v = pts[:, 1]-pts[:, 0], pts[:, 2]-pts[:, 0]
        total += .5*np.abs(u[:, 0]*v[:, 1]-u[:, 1]*v[:, 0]).sum()
    if abs(total-coverage.area) > remesh.area_tolerance(coverage):
        raise remesh.PatchError('Disk interior has inconsistent boundary coverage')
    return coverage, interfaces, cut_vertices


def grow_active(active, neighbours, seeds, rings):
    """Expand a local frontier without allocating a whole-domain index vector."""
    frontier = np.asarray(seeds, dtype=np.int64)
    added = []
    for _ in range(rings+1):
        frontier = frontier[~active[frontier]]
        if not len(frontier):
            break
        active[frontier] = True
        added.append(frontier)
        adjacent = neighbours[frontier].ravel()
        frontier = np.unique(adjacent[adjacent >= 0])
    return np.concatenate(added) if added else np.empty(0, dtype=np.int64)


def compact_region(vertices, triangles, active):
    # This allocation is proportional to the active seam, never the interiors.
    ids = np.concatenate([np.flatnonzero(active[lo:hi])+lo for lo, hi in batches(len(triangles))]) if len(triangles) else np.empty(0, dtype=np.int64)
    selected = triangles[ids]
    used, inverse = np.unique(selected, return_inverse=True)
    return vertices[used], inverse.reshape(-1, 3), ids


def boundary_index(vertices, triangles):
    if not len(triangles):
        return set(), {}
    neighbours = build_neighbours(triangles, len(vertices), vertices)
    faces, slots = np.nonzero(neighbours < 0)
    pairs = triangles[faces[:, None], (slots[:, None] + [1, 2]) % 3]
    points = {i: tuple(vertices[i, :2]) for i in np.unique(pairs)}
    edges = {tuple(sorted((points[a], points[b]))) for a, b in pairs}
    return edges, {xy: i for i, xy in points.items()}


def extend_cached_region(cached, inputs, pending):
    """Restore shared IDs along the previous disk cut when loading its next ring."""
    local_v, local_t, shared = cached
    vertices, triangles = [local_v], [local_t]
    keep = [np.ones(len(local_t), dtype=bool)]
    offset = len(local_v)
    for (v, t, _), (added, removed), cut_ids in zip(inputs, pending, shared):
        used, inverse = np.unique(t[added], return_inverse=True)
        ids = np.fromiter((cut_ids.get(i, -1) for i in used), dtype=np.int64, count=len(used))
        new = ids < 0
        new_vertices = v[used[new]]
        vertices.append(new_vertices)
        ids[new] = offset + np.arange(len(new_vertices))
        offset += len(new_vertices)
        triangles.append(ids[inverse].reshape(-1, 3))
        keep.append(~np.isin(added, removed))
    return np.vstack(vertices), np.vstack(triangles), np.concatenate(keep)


def stream_join(inputs, active, local_v, local_t, work):
    """Compact immutable interiors and append the patch using exact XY IDs.

    Input interiors occupy disjoint ownership boxes. Only patch coordinates
    can be shared across an interior/patch interface; their lookup is seam-sized.
    No global coordinate dictionary or global advanced-indexing copy is built.
    """
    mappings, counts = [], []
    for index, ((vertices, triangles, _), selection) in enumerate(zip(inputs, active)):
        ids = disk_array(work/f'ids-{index}.npy', (len(vertices),), np.int64)
        ids[:] = -1
        for lo, hi in batches(len(triangles)):
            ids[triangles[lo:hi][~selection[lo:hi]]] = 0
        count = 0
        for lo, hi in batches(len(ids)):
            used = np.flatnonzero(ids[lo:hi] >= 0)
            ids[lo+used] = sum(counts) + np.arange(count, count+len(used))
            count += len(used)
        counts.append(count)
        mappings.append(ids)
    lookup = {tuple(v[:2]): i for i, v in enumerate(local_v)}
    if len(lookup) != len(local_v):
        raise remesh.PatchError('Local merge contains coincident vertex IDs')
    local_ids = np.full(len(local_v), -1, dtype=np.int64)
    for (vertices, triangles, neighbours), selection, ids in zip(inputs, active, mappings):
        for _, endpoints, pair, _ in selected_boundary(vertices, triangles, neighbours, selection):
            for old, key in zip(endpoints, pair):
                local = lookup.get(key)
                if local is not None:
                    if local_ids[local] >= 0 and local_ids[local] != ids[old]:
                        raise remesh.PatchError('Immutable interiors share a vertex unexpectedly')
                    local_ids[local] = ids[old]
    added = np.flatnonzero(local_ids < 0)
    local_ids[added] = sum(counts) + np.arange(len(added))
    out_v = disk_array(work/'verts.npy', (sum(counts)+len(added), 3))
    nfaces = sum(sum(int(np.count_nonzero(~selection[lo:hi])) for lo, hi in batches(len(selection)))
                 for selection in active) + len(local_t)
    out_t = disk_array(work/'tris.npy', (nfaces, 3), np.int64)
    offset = 0
    for (vertices, triangles, _), selection, ids in zip(inputs, active, mappings):
        for lo, hi in batches(len(vertices)):
            used = lo+np.flatnonzero(ids[lo:hi] >= 0)
            out_v[ids[used]] = vertices[used]
        for lo, hi in batches(len(triangles)):
            block = ids[triangles[lo:hi][~selection[lo:hi]]]
            out_t[offset:offset+len(block)] = block
            offset += len(block)
    out_v[local_ids[added]] = local_v[added]
    for lo, hi in batches(len(local_t)):
        out_t[offset+lo:offset+hi] = local_ids[local_t[lo:hi]]
    out_v.flush(); out_t.flush()
    return out_v, out_t


def merge_disk_cavity(path_a, bounds_a, path_b, bounds_b, domain, constraints, options, output):
    started = time.perf_counter()
    # Keep failed native reproductions beside the final mesh even when staging its archive.
    directory = str(Path(options.get('out_prefix') or output).parent)
    width = options.get('remesh_band_width')
    width = float(2*math.sqrt(options['max_area']) if width is None else width)
    limit = int(options.get('remesh_max_attempts', 12))
    max_active = int(options.get('remesh_max_active_faces', 2000000))
    if limit < 1 or max_active < 1 or not np.isfinite(width) or width <= 0:
        raise remesh.PatchError('Seam width, attempt count and active-face limit must be positive')
    axis, seam = remesh._shared_interface(bounds_a, bounds_b)
    bounds = [min(bounds_a[0], bounds_b[0]), min(bounds_a[1], bounds_b[1]),
              max(bounds_a[2], bounds_b[2]), max(bounds_a[3], bounds_b[3])]
    target = remesh._clipped_domain(domain, bounds, 0, width/2, options.get('out_prefix'),
                                    options.get('mesher_debug', False))
    with tempfile.TemporaryDirectory(prefix='seam-store-', dir=directory) as name, ExitStack() as stack:
        work = Path(name)
        inputs, active = [], []
        input_seconds = 0.0
        for index, (path, core) in enumerate(((path_a, bounds_a), (path_b, bounds_b))):
            read_started = time.perf_counter()
            v, t = stack.enter_context(open_mesh_arrays(path, work))
            input_seconds += time.perf_counter() - read_started
            neighbours = (build_neighbours(t, len(v), v, chunk_size=32768,
                output_path=work/f'neighbours-{index}.npy', temporary_directory=work)
                if len(t) else np.empty((0, 3), dtype=np.int64))
            selection = disk_array(work/f'active-{index}.npy', (len(t),), bool)
            for lo, hi in batches(len(t)):
                selection[lo:hi] = ~remesh._retain_mask(v, t[lo:hi], core, axis, seam, width)
            # A four-ring collar gives local parent repairs room before touching disk interiors.
            seeds = []
            for lo, hi in batches(len(t)):
                adjacent = neighbours[lo:hi]
                touches_active = ((adjacent >= 0) & selection[adjacent]).any(axis=1)
                seeds.extend(lo + np.flatnonzero(~selection[lo:hi] & touches_active))
            grow_active(selection, neighbours, seeds, 3)
            inputs.append((v, t, neighbours)); active.append(selection)
        reports = []
        cached = None
        pending = None
        total = sum(len(t) for _, t, _ in inputs)
        peak_active = 0
        reused = 0
        for attempt in range(limit):
            coverage, cuts, cut_vertices = [], [], []
            for values, selection in zip(inputs, active):
                polygon, interfaces, endpoints = inactive_coverage(*values, selection)
                coverage.append(polygon)
                cuts.append(interfaces)
                cut_vertices.append(endpoints)
            outside = remesh.polygon_domain(union_all(coverage))
            if outside.difference(target).area > remesh.area_tolerance(target):
                raise remesh.PatchError('Immutable interiors extend beyond the merge domain')
            local_domain = remesh.polygon_domain(target.difference(outside))
            selected_count = sum(int(np.count_nonzero(selection)) for selection in active)
            active_count = (selected_count if cached is None
                            else len(cached[1]) + sum(len(added) for added, _ in pending))
            peak_active = max(peak_active, active_count)
            if active_count > max_active:
                raise remesh.PatchError(f'Seam active-face budget exceeded: {active_count} > {max_active}; '
                                       'increase mpi_remesh_max_active_faces only if memory permits')
            if cached is None:
                av, at, _ = compact_region(*inputs[0][:2], active[0])
                bv, bt, _ = compact_region(*inputs[1][:2], active[1])
                initial_keep = None
            else:
                reused += len(cached[1])
                av, at, initial_keep = extend_cached_region(cached, inputs, pending)
                # The expanded seam is already indexed as a single working mesh.
                bv, bt = av[:0], at[:0]
            debug_log(options, f'Seam region [{Path(output).stem}]: active={active_count}, '
                      f'interiors_on_disk={total-selected_count}')
            try:
                lv, lt, report = remesh.merge_cavity(av, at, bounds_a, bv, bt, bounds_b,
                    domain, constraints, options, directory, prepared_domain=local_domain, initial_keep=initial_keep)
                # Local checks preserve the target footprint; verify the disk-facing
                # interface by exact endpoints before treating it as a completed merge.
                if len(lt):
                    local_coverage = remesh.mesh_coverage(lv, lt)
                    if local_coverage.symmetric_difference(local_domain).area > remesh.area_tolerance(local_domain):
                        raise remesh.PatchError('Seam region has inconsistent target coverage')
                elif not local_domain.is_empty:
                    raise remesh.PatchError('Seam region is unexpectedly empty')
                keys, boundary_ids = boundary_index(lv, lt)
                affected = [sorted({face for edge, face in interfaces.items() if edge not in keys})
                            for interfaces in cuts]
                reports.append(report)
                if any(affected):
                    shared = [{i: boundary_ids[xy] for i, xy in endpoints.items() if xy in boundary_ids}
                              for endpoints in cut_vertices]
                    cached = lv, lt, shared
            except remesh.PatchError as exc:
                cause = exc
                while cause is not None:
                    if (isinstance(cause, remesh.BackendTerminated) or 'Patch vertex budget exceeded' in str(cause)
                            or 'Seam active-face budget exceeded' in str(cause)):
                        raise
                    cause = cause.__cause__
                # A strict quality failure at an artificial outer cut can require
                # more neighbours. Never relax the criteria to keep an interior.
                affected = [sorted(set(interfaces.values())) for interfaces in cuts]
                if not any(affected) or attempt+1 == limit:
                    raise
                cached = None
                debug_log(options, f'Expanding seam region after rejected local repair: {exc}')
            else:
                if not any(affected):
                    out_v, out_t = stream_join(inputs, active, lv, lt, work)
                    if len(out_t):
                        neighbours = build_neighbours(out_t, len(out_v), out_v, 32768,
                            work/'final-neighbours.npy', work)
                        check_boundary_geometry(out_v, out_t, neighbours, work)
                    # Copy NPY members directly into the backward-compatible NPZ archive.
                    write_started = time.perf_counter()
                    temporary = work/'complete.npz'
                    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_STORED, allowZip64=True) as archive:
                        archive.write(work/'verts.npy', 'verts.npy')
                        archive.write(work/'tris.npy', 'tris.npy')
                    os.replace(temporary, output)
                    report.update(seam_attempts=attempt+1, seam_attempt_reports=reports[:-1],
                        peak_loaded_faces=peak_active, seam_work_reused=reused,
                        boundary_halo_area=target.area-domain.intersection(box(*bounds)).area,
                        disk_interior_faces=len(out_t)-len(lt), output_faces=len(out_t),
                        input_seconds=input_seconds, output_npz_seconds=time.perf_counter()-write_started,
                        elapsed_seconds=time.perf_counter()-started)
                    return report
            if attempt+1 == limit:
                raise remesh.PatchError('Seam-region expansion exhausted before the disk interface matched')
            pending = []
            for selection, (_, _, neighbours), seeds in zip(active, inputs, affected):
                # Remove the failed interface parents and their immediate ring;
                # load three further rings as retained repair candidates.
                if seeds:
                    adjacent = neighbours[seeds].ravel()
                    removed = np.unique(np.concatenate((seeds, adjacent[adjacent >= 0])))
                    removed = removed[~selection[removed]]
                    added = grow_active(selection, neighbours, seeds, 4)
                else:
                    removed = added = np.empty(0, dtype=np.int64)
                pending.append((added, removed))
            if reports:
                reports[-1]['active_faces'] = active_count
        raise remesh.PatchError('No conforming seam-region replacement')
