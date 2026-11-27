"""Split a fully released seam without relaxing CGAL's physical or quality constraints.

Each artificial cut has one ordered subdivision shared by both incident pieces.
A worker may discover extra vertices there. Retry just the side missing a point,
then remove a troublesome cut locally if the two triangulations keep disagreeing.
Only validated pieces reach the restart cache. The pair's ordinary disk-interface
check still decides whether the assembled candidate can be published.
"""
from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import zipfile
import zlib

import numpy as np
from shapely import STRtree
from shapely.geometry import LineString, MultiPolygon, box
from shapely.errors import GEOSException

from mesher.mesh_topology import build_neighbours, MeshTopologyError
from . import mpi_remesh as remesh
from .mpi_progress import report_phase


@dataclass
class Piece:
    # Reclip the original polygon when removing a cut. Unioning clipped halves
    # would leave an artificial vertex on a formerly straight physical boundary.
    source: object
    lo: float
    hi: float
    domain: object
    result: object = None


def _slab(source, axis, lo, hi):
    bounds = list(source.bounds)
    if lo == bounds[axis] and hi == bounds[axis+2]:
        return Piece(source, lo, hi, source)
    bounds[axis], bounds[axis+2] = lo, hi
    clipped = source.intersection(box(*bounds))
    if clipped.area > 0:
        return Piece(source, lo, hi, remesh.polygon_domain(clipped))
    return None


def _pieces(domain, axis, count):
    # Approximate work using area within the active region. Keep holes with
    # their polygon and batch tiny islands so each CGAL start has useful work.
    target = domain.area / count
    small, area, pieces = [], 0., []

    def flush():
        if small:
            source = small[0] if len(small) == 1 else MultiPolygon(small)
            pieces.append(_slab(source, axis, source.bounds[axis], source.bounds[axis+2]))
            small.clear()

    for polygon in getattr(domain, 'geoms', [domain]):
        n = math.ceil(polygon.area / target)
        if n == 1:
            if area + polygon.area > target:
                flush()
                area = 0.
            small.append(polygon)
            area += polygon.area
            continue
        lo, hi = polygon.bounds[axis], polygon.bounds[axis+2]
        cuts = np.linspace(lo, hi, n+1)
        for a, b in zip(cuts, cuts[1:]):
            piece = _slab(polygon, axis, a, b)
            if piece is not None:
                pieces.append(piece)
    flush()
    return pieces


def _cuts(pieces):
    # Physical boundaries can run along the same coordinate as a cut (holes,
    # inlets). Only positive-length boundary shared by the two slabs belongs
    # to this negotiation. Point contacts share a vertex during exact assembly.
    cuts = {}
    for i, (left, right) in enumerate(zip(pieces, pieces[1:])):
        if left.source is right.source and left.hi == right.lo:
            shared = left.domain.boundary.intersection(right.domain.boundary)
            for line in remesh.line_parts(shared):
                xy = list(line.coords)
                for a, b in zip(xy, xy[1:]):
                    if a != b:
                        cuts[tuple(sorted((a, b)))] = (i, i+1)
    return cuts


def _solve_piece(job, options):
    """Publish a complete, independently validated native result in one rename."""
    domain, streams, pins, cache = job
    # checkpoint_run and parent generations scope the directory. The key also
    # changes when the exact footprint, clipped streams or requested pins change.
    digest = hashlib.sha256(domain.normalize().wkb)
    for line in streams:
        digest.update(line.wkb)
    digest.update(repr(pins).encode())
    path = cache / (digest.hexdigest() + '.npz')
    if path.exists():
        try:
            with np.load(path, allow_pickle=False) as data:
                # Reading the payload checks its ZIP CRC before trusting a
                # previous allocation's result. No incomplete writer is visible.
                for key in ('vertices', 'triangles', 'boundary', 'constrained'):
                    data[key]  # np.load is lazy, access each array to check its CRC.
                chains = json.loads(data['chains'].item())
            return path, chains, {}, True
        except (OSError, ValueError, EOFError, KeyError, zipfile.BadZipFile, zlib.error):
            pass

    segments = [LineString([a, b]) for chain in pins for a, b in zip(chain, chain[1:])]
    directory = Path(options['out_prefix']).parent if 'out_prefix' in options else cache.parent
    patch = remesh.remesh_domain(domain, streams, options, directory, segments, return_constrained_edges=True)
    vertices, triangles = patch
    # Confirm that every reported child is an actual boundary edge. Merely
    # having the same vertex set would allow a hanging node or an interior edge.
    adjacency = build_neighbours(triangles, len(vertices), vertices)
    faces, slots = np.nonzero(adjacency < 0)
    boundary = set(zip(triangles[faces, (slots+1)%3], triangles[faces, (slots+2)%3]))
    lineage = {item['id']: item['vertices'] for item in patch.interfaces}
    if set(lineage) != set(range(len(segments))):
        raise remesh.PatchError('Missing cavity cut lineage')
    chains = []
    index = 0
    for chain in pins:
        agreed = [chain[0]]
        for a, b in zip(chain, chain[1:]):
            ids = lineage[index]
            index += 1
            if len(ids) < 2 or min(ids) < 0 or max(ids) >= len(vertices):
                raise remesh.PatchError('Invalid cavity cut vertex IDs')
            xy = list(map(tuple, vertices[ids, :2].tolist()))
            axis = 0 if a[0] == b[0] else 1
            if (xy[0] != a or xy[-1] != b or xy != sorted(set(xy)) or any(p[axis] != a[axis] for p in xy)
                    or any((j, k) not in boundary and (k, j) not in boundary for j, k in zip(ids, ids[1:]))):
                raise remesh.PatchError('Invalid cavity cut subdivision')
            agreed.extend(xy[1:])
        chains.append(agreed)
    with tempfile.NamedTemporaryFile(dir=cache, suffix='.npz', delete=False) as stream:
        tmp = Path(stream.name)
        try:
            border = np.unique(triangles[faces[:, None], (slots[:, None] + [1, 2]) % 3])
            constrained = np.asarray(patch.constrained_edges or [], dtype=np.int64).reshape(-1, 2)
            np.savez(stream, vertices=vertices, triangles=triangles, chains=json.dumps(chains),
                     boundary=border, constrained=constrained)
            stream.close()
            os.replace(tmp, path)
        finally:
            tmp.unlink(missing_ok=True)
    return path, chains, patch.stats, False


class CavityTeam:
    """A flat MPI team for one tile pair. Rank zero also runs patch jobs."""
    def __init__(self, comm, options):
        self.comm = comm
        self.options = options

    def _work(self, jobs):
        results = []
        for job in jobs:
            try:
                result = _solve_piece(job, self.options)
                results.append(('ok', result))
            except remesh.BackendTerminated as exc:
                results.append(('fatal', str(exc)))
            except (MeshTopologyError, GEOSException) as exc:
                results.append(('geometry', str(exc)))
            except Exception as exc:
                # All ranks must reach gather even on I/O or backend failures.
                results.append(('error', f'{type(exc).__name__}: {exc}'))
        return results

    def serve(self):
        while True:
            jobs = self.comm.scatter(None, root=0)
            if jobs is None:
                return
            self.comm.gather(self._work(jobs), root=0)

    def stop(self):
        self.comm.scatter([None] * self.comm.size, root=0)

    def map(self, jobs):
        # Costlier regions start first. The on-disk outputs keep triangle
        # buffers off the MPI wire and out of the coordinator's result queue.
        order = sorted(range(len(jobs)), key=lambda i: jobs[i][0].area, reverse=True)
        buckets = [order[i::self.comm.size] for i in range(self.comm.size)]
        local = self.comm.scatter([[jobs[i] for i in bucket] for bucket in buckets], root=0)
        gathered = self.comm.gather(self._work(local), root=0)
        results = [None] * len(jobs)
        for bucket, batch in zip(buckets, gathered):
            for i, (status, value) in zip(bucket, batch):
                if status == 'fatal':
                    raise remesh.BackendTerminated(value)
                if status == 'error':
                    raise RuntimeError(value)
                results[i] = None if status == 'geometry' else value
        return results


def _assemble(pieces, domain, constraints, options):
    # Read each archive once into the combined working arrays. Deduplicate by
    # exact XY in one pass rather than repeatedly copying a growing mesh.
    vertices, triangles, constrained = [], [], []
    shared_ids = {}
    offset = 0
    for piece in pieces:
        with np.load(piece.result[0]) as data:
            v, t = data['vertices'], data['triangles']
            edge_ids, border = data['constrained'], data['boundary']
        # Only boundary vertices can belong to two nonoverlapping pieces.
        # Keep this map proportional to perimeter and enforce the vertex budget
        # before accumulating another piece in the coordinator.
        ids = np.full(len(v), -1, dtype=np.int64)
        for j in border:
            ids[j] = shared_ids.get(tuple(v[j, :2]), -1)
        new = ids < 0
        count = int(new.sum())
        if offset + count > options.get('remesh_max_vertices', 500000):
            raise remesh.PatchError('Patch vertex budget exceeded across cavity pieces')
        ids[new] = offset + np.arange(count)
        for j in border[new[border]]:
            shared_ids[tuple(v[j, :2])] = ids[j]
        vertices.append(v[new])
        triangles.append(ids[t])
        constrained.append(ids[edge_ids])
        offset += count
    v, t = np.concatenate(vertices), np.concatenate(triangles)
    neighbours = build_neighbours(t, len(v), v)
    coverage = remesh.mesh_coverage(v, t, neighbours)
    if coverage.symmetric_difference(domain).area > remesh.area_tolerance(domain):
        raise remesh.PatchError('Parallel cavity assembly changed the domain')
    if constraints:
        remesh._validate_stream_edges(v, t, neighbours, domain, constraints, np.concatenate(constrained))
    return remesh.PatchResult(v, t, [], {})


def solve(domain, constraints, options, axis, face_count, phase):
    """Mesh one released cavity using its face count as a cost estimate.

    The axis runs along the original tile interface. Each key in cuts is an
    oriented segment, always ordered by increasing coordinate, and its values
    identify the two incident pieces. Pins contain the ordered vertices each
    future solve must reproduce on that segment. A result stays eligible as
    long as it already supplies all the pins requested by its neighbours.
    """
    count = math.ceil(face_count / options['cavity_target_faces'])
    try:
        pieces = _pieces(domain, axis, count)
    except (MeshTopologyError, GEOSException) as exc:
        remesh.debug_log(options, f'Cannot partition cavity: {exc}')
        return None
    team = options['_cavity_team']
    cache = options['_cavity_cache']
    cache.mkdir(parents=True, exist_ok=True)
    tree = STRtree(constraints)
    cuts = _cuts(pieces)
    pins = {edge: edge for edge in cuts}
    rounds = defaultdict(int)
    stats = defaultdict(float)
    hits = 0
    solves = 0
    while len(pieces) > 1:
        edges = [[] for _ in pieces]
        for edge, owners in cuts.items():
            for i in owners:
                edges[i].append(edge)
        jobs, indices = [], []
        for i, piece in enumerate(pieces):
            if piece.result is not None:
                continue
            streams = [part for k in sorted(tree.query(piece.domain))
                       for part in remesh.line_parts(constraints[k].intersection(piece.domain))]
            jobs.append((piece.domain, streams, [pins[e] for e in edges[i]], cache))
            indices.append(i)
        report_phase(options, f'{phase} | solving {len(jobs)} of {len(pieces)} cavity pieces '
                             f'on {team.comm.size} ranks')
        results = team.map(jobs)
        remove = set()
        observed = {}
        for i, result in zip(indices, results):
            if result is None:
                # A cut can introduce an infeasible acute angle. Restore the
                # original local polygon before assigning blame to a stream.
                for j in (i-1, i):
                    if 0 <= j < len(pieces)-1 and pieces[j].source is pieces[j+1].source:
                        remove.add(j)
                if not remove.intersection({i-1, i}):
                    return None
            else:
                path, chains, timings, cached = result
                hits += int(cached)
                solves += int(not cached)
                for key, value in timings.items():
                    if key.endswith('_seconds') or key in ('refinement_steps', 'quality_repair_vertices', 'raster_evaluations'):
                        stats[key] += value
                pieces[i].result = path, {e: tuple(map(tuple, ch)) for e, ch in zip(edges[i], chains)}
        if not remove:
            for edge, (i, j) in cuts.items():
                a, b = pieces[i].result[1][edge], pieces[j].result[1][edge]
                if a == b:
                    continue
                rounds[edge] += 1
                # Allow two rounds of incorporating new pins after the first
                # disagreement. A persistent split cascade is cheaper to solve
                # by removing this artificial cut and keeping the native quality checks.
                if rounds[edge] >= 3:
                    remove.add(i)
                else:
                    wanted = tuple(sorted(set(a) | set(b)))
                    pins[edge] = wanted
                    observed[i] = observed.get(i, False) or a != wanted
                    observed[j] = observed.get(j, False) or b != wanted
            if not remove:
                if not observed:
                    try:
                        result = _assemble(pieces, domain, constraints, options)
                    except (MeshTopologyError, GEOSException) as exc:
                        if 'Patch vertex budget exceeded' in str(exc):
                            raise
                        remesh.debug_log(options, f'Parallel cavity join rejected: {exc}')
                        return None
                    result.stats.update(stats, parallel_pieces=len(pieces), cavity_ranks=team.comm.size, cached_piece_reads=hits,
                                        piece_solves=solves)
                    return result
                for i, retry in observed.items():
                    if retry:
                        pieces[i].result = None
                continue
        remesh.debug_log(options, f'Removing {len(remove)} artificial cavity cuts after local refinement')
        # Only neighbours of a failed interface lose their accepted piece.
        # Recompute incidence after removing cuts but keep all other archives.
        next_pieces = []
        i = 0
        while i < len(pieces):
            j = i
            while j in remove:
                j += 1
            next_pieces.append(pieces[i] if i == j else
                               _slab(pieces[i].source, axis, pieces[i].lo, pieces[j].hi))
            i = j + 1
        pieces = next_pieces
        cuts = _cuts(pieces)
        pins = {e: pins.get(e, e) for e in cuts}
        for i, piece in enumerate(pieces):
            if piece.result is not None and set(piece.result[1]) != {e for e, owners in cuts.items() if i in owners}:
                piece.result = None
    return None
