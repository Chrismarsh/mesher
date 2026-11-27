"""Run identity and bounded inspection of committed MPI tile/merge archives."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import uuid
import zipfile
import zlib

import numpy as np
from osgeo import gdal


# Hold one exclusive lock across cleanup and publication so two jobs cannot mix their tile
# generations in the same output directory.
@contextmanager
def output_lock(directory):
    import fcntl
    directory = Path(directory)
    directory.parent.mkdir(parents=True, exist_ok=True)
    # A sibling survives cleanup of a fresh run. Never unlink a live lock inode.
    with (directory.parent / f'.{directory.name}.mesher.lock').open('a') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Another Mesher run is using {directory}') from exc
        yield


def input_files(paths):
    """Identify local datasets, including shapefile members and nested VRT sources."""
    files = {}
    # Key by real path so files reached through several VRTs are visited once. Returning
    # False also stops traversal from following a reference back to a known dataset.
    def record(path):
        key = os.path.realpath(path)
        if key in files:
            return False
        try:
            stat = os.stat(path)
        except OSError:
            files[key] = None
        else:
            files[key] = [stat.st_size, stat.st_mtime_ns]
        return True

    pending = list(map(os.fspath, paths))
    while pending:
        path = pending.pop()
        if not record(path):
            continue
        dataset = gdal.OpenEx(path)
        if dataset is None:
            raise ValueError(f'Cannot open meshing input {path}')
        for member in dataset.GetFileList() or ():
            if member.lower().endswith('.vrt'):
                pending.append(member)
            else:
                record(member)
    return files


# A normal resume needs a validated checkpoint. Legacy mode also accepts an absent plan so the
# caller can inspect older tile markers.
def read_plan(path, mode):
    if not mode:
        return None
    if not path.exists():
        if mode == 'legacy':
            return None
        raise ValueError(f'No MPI checkpoint at {path}. For an older interrupted run, '
                         'use --resume-legacy with the original configuration and inputs.')
    with path.open() as stream:
        plan = json.load(stream)
    if (not isinstance(plan, dict) or plan.get('version') != 1 or not plan.get('run_id')
            or not isinstance(plan.get('recipe'), dict) or not {'sources', 'parameters'} <= plan['recipe'].keys()
            or type(plan.get('legacy')) is not bool or not isinstance(plan.get('grid'), list)
            or len(plan['grid']) != 2 or any(type(n) is not int or n < 1 for n in plan['grid'])):
        raise ValueError(f'Invalid MPI checkpoint: {path}')
    return plan


# Write beside the destination and rename the complete JSON into place so readers cannot observe
# a partial checkpoint.
def write_marker(path, data):
    temporary = path.with_name(f'.{path.name}.{uuid.uuid4().hex}.tmp')
    try:
        with temporary.open('x') as stream:
            json.dump(data, stream)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


# A restart may reuse a run ID only when its recipe and tile grid match. A fresh run gets a new
# generation of checkpoints.
def start_plan(path, previous, recipe, grid, legacy):
    recipe = json.loads(json.dumps(recipe))  # normalize Python tuples for comparison on reload
    if previous:
        differences = [key for key in recipe if previous['recipe'].get(key) != recipe[key]]
        if differences or previous['grid'] != list(grid):
            raise ValueError('MPI checkpoint is incompatible with the requested ' +
                             ', '.join(differences or ['tile grid']))
        return previous
    plan = dict(version=1, run_id=uuid.uuid4().hex, grid=list(grid), recipe=recipe, legacy=legacy)
    write_marker(path, plan)
    return plan


def archive_info(path):
    """Read ZIP directory and NPY headers only, including for legacy archives."""
    members = []
    with zipfile.ZipFile(path) as archive:
        if sorted(archive.namelist()) != ['tris.npy', 'verts.npy']:
            raise ValueError('Expected verts.npy and tris.npy')
        for name, kinds in [('verts.npy', 'f'), ('tris.npy', 'iu')]:
            info = archive.getinfo(name)
            # A large restart may inspect thousands of archives. The header
            # and ZIP directory let us check the array dimensions and compare
            # the recorded size and CRC without scanning the triangle payload.
            with archive.open(info) as stream:
                version = np.lib.format.read_magic(stream)
                if version == (1, 0):
                    shape, _, dtype = np.lib.format.read_array_header_1_0(stream)
                elif version == (2, 0):
                    shape, _, dtype = np.lib.format.read_array_header_2_0(stream)
                else:
                    raise ValueError(f'Unsupported NPY version {version}')
                if (len(shape) != 2 or shape[1] != 3 or dtype.kind not in kinds
                        or info.file_size != stream.tell() + shape[0]*3*dtype.itemsize):
                    raise ValueError(f'Invalid array {name} in {path}')
            members.append([info.file_size, info.CRC])
    return members, shape[0]  # tris is the second member


# Read both inputs in task order. The child stores their generation IDs in this same order when
# recording which pair it merged.
def parent_markers(task):
    for tile in (task['tile_a'], task['tile_b']):
        with open(tile['meta']) as stream:
            yield json.load(stream)


def publish(staged, prefix, marker, task):
    """The JSON rename commits the archive. A killed writer leaves a missing marker."""
    archive = Path(str(prefix) + '.npz')
    metadata = Path(str(prefix) + '.json')
    if 'checkpoint_run' in task:
        members, _ = archive_info(staged)
        generations = [m['checkpoint']['generation'] for m in parent_markers(task)] if 'tile_a' in task else []
        marker['checkpoint'] = dict(run=task['checkpoint_run'], generation=uuid.uuid4().hex,
                                    archive=members, parents=generations)
    staged_marker = Path(staged).with_suffix('.json')
    with staged_marker.open('w') as stream:
        json.dump(marker, stream)
    # Retire the old marker first. If this worker dies between renames, a
    # restart must schedule the pair again even if its NPZ already exists.
    metadata.unlink(missing_ok=True)
    os.replace(staged, archive)
    os.replace(staged_marker, metadata)


# Schedule any tile or merge whose committed archive or parent generations no longer match. A
# rebuilt parent must invalidate previously completed descendants.
def unfinished(tasks, plan, label):
    pending = []
    for task in tasks:
        merge = 'tile_a' in task
        prefix = task['out_prefix'] if merge else task['base_dir'] + task['tile_prefix']
        try:
            with open(prefix + '.json') as stream:
                marker = json.load(stream)
            if merge:
                records = list(parent_markers(task))
                a, b = (parent['core_bbox'] for parent in records)
                bounds = [min(a[0], b[0]), min(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3])]
            else:
                bounds = task['tile_bbox']
            if (marker['merge_method'] != 'remesh' or marker['core_bbox'] != list(bounds)
                    or marker['band_width'] != (0 if merge else task['band_width'])):
                raise ValueError('Different tile or merge geometry')
            members, faces = archive_info(prefix + '.npz')
            if marker['empty_tile'] != (faces == 0) or (merge and marker['remesh_report']['output_faces'] != faces):
                raise ValueError('Face count mismatch')
            # Rebuilding a parent can change connectivity while preserving its
            # filename, bounds and face count. Its new generation invalidates
            # the merges that used it, leaving sibling branches reusable.
            generations = [m['checkpoint']['generation'] for m in records] if merge else []
            checkpoint = marker.get('checkpoint')
            if checkpoint is None:
                valid = plan['legacy'] and all(key.startswith('legacy:') for key in generations)
            else:
                valid = (checkpoint['run'] == plan['run_id'] and checkpoint['archive'] == members
                         and checkpoint['parents'] == generations and isinstance(checkpoint['generation'], str)
                         and bool(checkpoint['generation']))
        except (OSError, EOFError, ValueError, KeyError, TypeError, AttributeError, zlib.error, zipfile.BadZipFile):
            valid = False
        if not valid:
            pending.append(task)
        elif checkpoint is None:
            marker['checkpoint'] = dict(run=plan['run_id'], generation='legacy:' + uuid.uuid4().hex,
                                        archive=members, parents=generations)
            write_marker(Path(prefix + '.json'), marker)
    print(f'{label}: reusing {len(tasks)-len(pending)}/{len(tasks)}, scheduling {len(pending)}', flush=True)
    return pending


def check_legacy_grid(directory, base, tasks):
    """Older markers identify their tile and band, but cannot authenticate settings."""
    found = False
    prefix = base + '_mpi_tile_'
    for path in Path(directory).glob(prefix + '*.json'):
        tile_id = path.stem[len(prefix):]
        if not tile_id.isdecimal():
            continue
        try:
            marker = json.loads(path.read_text())
        except (ValueError, OSError):
            continue
        if not isinstance(marker, dict) or 'core_bbox' not in marker:
            continue
        index = int(tile_id)
        if (marker.get('merge_method') != 'remesh' or marker.get('checkpoint') is not None
                or index >= len(tasks) or marker['core_bbox'] != list(tasks[index]['tile_bbox'])
                or marker.get('band_width') != tasks[index]['band_width']):
            raise ValueError(f'{path} does not match the original legacy tile grid or remesh band. '
                             'Use the original layout and configuration.')
        found = True
    if not found:
        raise ValueError('No identifiable legacy MPI tile markers in ' + str(directory))
    print('Adopting legacy MPI outputs using the supplied original inputs and tolerances.', flush=True)
