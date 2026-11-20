"""Disk-backed mesh IO. Only bounded numeric or JSON batches enter the heap."""

from contextlib import closing, contextmanager
import json
import os
from pathlib import Path
import shutil
import sqlite3
import struct
import tempfile
from xml.sax.saxutils import quoteattr
import zipfile

import numpy as np

from .mesh_topology import MeshTopologyError

BATCH_SIZE = 32768


def batches(count, size=BATCH_SIZE):
    if size < 1:
        raise ValueError('Batch size must be positive')
    for start in range(0, count, size):
        yield start, min(count, start + size)


def disk_array(path, shape, dtype=float):
    return np.lib.format.open_memmap(path, mode='w+', dtype=dtype, shape=shape)


def copy_array(path, source, dtype=None):
    target = disk_array(path, np.shape(source), dtype or np.asarray(source).dtype)
    for lo, hi in batches(len(source)):
        target[lo:hi] = source[lo:hi]
    target.flush()
    return target


@contextmanager
def open_mesh_arrays(path, directory=None):
    """Unpack old and current NPZ files incrementally, then map their NPY members.

    np.load(..., mmap_mode=...) does not map NPZ members. Explicit extraction
    also supports compressed legacy archives without materializing an array.
    """
    with tempfile.TemporaryDirectory(prefix='mesh-arrays-', dir=directory) as work:
        arrays = {}
        with zipfile.ZipFile(path) as archive:
            for key in ('verts', 'tris'):
                target = Path(work) / (key + '.npy')
                with archive.open(key + '.npy') as source, target.open('wb') as out:
                    shutil.copyfileobj(source, out, length=1024 * 1024)
                arrays[key] = np.load(target, mmap_mode='r', allow_pickle=False)
        try:
            yield arrays['verts'], arrays['tris']
        finally:
            for array in arrays.values():
                array._mmap.close()


def boundary_edges(triangles, neighbours):
    for lo, hi in batches(len(triangles)):
        faces, slots = np.nonzero(neighbours[lo:hi] < 0)
        faces += lo
        for face, slot in zip(faces, slots):
            yield int(face), int(slot), int(triangles[face, (slot+1) % 3]), int(triangles[face, (slot+2) % 3])


def check_boundary_geometry(vertices, triangles, neighbours, directory, tolerance=1e-8):
    """Check the same boundary coincidences/T-junctions using an on-disk index.

    SQLite R-tree bounds are conservative float32 boxes; actual tests always use
    the stored float64 coordinates. No snapping or tolerance-based welding.
    """
    with tempfile.TemporaryDirectory(prefix='boundary-index-', dir=directory) as work:
        with closing(sqlite3.connect(str(Path(work) / 'index.sqlite'))) as db:
            db.execute('PRAGMA cache_size=-8192')
            db.execute('PRAGMA temp_store=FILE')
            db.execute('CREATE TABLE points (id INTEGER PRIMARY KEY, x REAL, y REAL)')
            db.execute('CREATE VIRTUAL TABLE spatial USING rtree(id, x0, x1, y0, y1)')
            for _, _, a, b in boundary_edges(triangles, neighbours):
                for vid in (a, b):
                    x, y = map(float, vertices[vid, :2])
                    if db.execute('INSERT OR IGNORE INTO points VALUES (?,?,?)', (vid, x, y)).rowcount:
                        db.execute('INSERT INTO spatial VALUES (?,?,?,?,?)', (vid, x, x, y, y))
            db.commit()

            def nearby(x0, x1, y0, y1):
                return db.execute('SELECT p.id,p.x,p.y FROM spatial s JOIN points p ON p.id=s.id '
                                  'WHERE s.x0<=? AND s.x1>=? AND s.y0<=? AND s.y1>=?',
                                  (x1, x0, y1, y0))

            for vid, x, y in db.execute('SELECT * FROM points'):
                for other, xx, yy in nearby(x-tolerance, x+tolerance, y-tolerance, y+tolerance):
                    if other != vid and np.hypot(xx-x, yy-y) <= tolerance:
                        raise MeshTopologyError('Coincident boundary vertices in MPI mesh')
            for face, _, a, b in boundary_edges(triangles, neighbours):
                p, q = vertices[a, :2], vertices[b, :2]
                dx, dy = q-p
                length = np.hypot(dx, dy)
                for vid, x, y in nearby(min(p[0], q[0])-tolerance, max(p[0], q[0])+tolerance,
                                       min(p[1], q[1])-tolerance, max(p[1], q[1])+tolerance):
                    if vid in (a, b):
                        continue
                    along = ((x-p[0])*dx + (y-p[1])*dy) / length
                    distance = abs(dx*(y-p[1])-dy*(x-p[0])) / length
                    if tolerance < along < length-tolerance and distance <= tolerance:
                        raise MeshTopologyError(f'Hanging boundary nodes at face {face}')


def dump_json(value, stream):
    """Write the existing CHM JSON schema without whole-array list conversion."""
    if isinstance(value, dict):
        stream.write('{')
        for index, (key, item) in enumerate(value.items()):
            if index:
                stream.write(',')
            stream.write(json.dumps(key) + ':')
            dump_json(item, stream)
        stream.write('}')
    elif isinstance(value, np.ndarray):
        stream.write('[')
        first = True
        for lo, hi in batches(len(value)):
            for row in value[lo:hi].tolist():
                if not first:
                    stream.write(',')
                stream.write(json.dumps(row, separators=(',', ':')))
                first = False
        stream.write(']')
    else:
        json.dump(value, stream)


def write_vtu(path, mesh, parameters, initial_conditions):
    """Stream one VTU piece with appended binary arrays and global point IDs."""
    m = mesh['mesh']
    vertices, triangles = m['vertex'], m['elem']
    n, nv = len(triangles), len(vertices)
    arrays, offset = [], 0

    def descriptor(name, dtype, count, components, blocks):
        nonlocal offset
        dt = np.dtype(dtype)
        vtk_type = {('f', 8): 'Float64', ('i', 8): 'Int64', ('u', 1): 'UInt8'}[dt.kind, dt.itemsize]
        tag = (f'<DataArray type="{vtk_type}" Name={quoteattr(name)} NumberOfComponents="{components}" '
               f'format="appended" offset="{offset}"/>')
        size = count * components * dt.itemsize
        arrays.append((dt, size, blocks))
        offset += 8 + size
        return tag

    def slices(a):
        return lambda: (a[lo:hi] for lo, hi in batches(len(a)))

    points = descriptor('Points', '<f8', nv, 3, slices(vertices))
    cells = descriptor('connectivity', '<i8', n*3, 1, slices(triangles))
    cells += descriptor('offsets', '<i8', n, 1,
                        lambda: (np.arange(lo+1, hi+1, dtype=np.int64)*3 for lo, hi in batches(n)))
    cells += descriptor('types', 'u1', n, 1, lambda: (np.full(hi-lo, 5, dtype=np.uint8) for lo, hi in batches(n)))
    data = descriptor('elevation', '<f8', n, 1,
                      lambda: (np.asarray(vertices)[np.asarray(triangles[lo:hi]), 2].mean(axis=1) for lo, hi in batches(n)))
    data += descriptor('cellid', '<i8', n, 1,
                       lambda: (np.arange(lo, hi, dtype=np.int64) for lo, hi in batches(n)))
    data += descriptor('area', '<f8', n, 1, lambda: (np.ones(hi-lo) for lo, hi in batches(n)))
    for prefix, values in (('[param] ', parameters), ('[ic] ', initial_conditions)):
        for key, array in values.items():
            if len(array) != n:
                raise ValueError(f'{key}: expected {n} cell values, got {len(array)}')
            data += descriptor(prefix + key, '<f8', n, 1, slices(array))
    projection = ' '.join(str(c) for c in m.get('proj4', '').encode('utf-8') + b'\0')
    header = ('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="1.0" '
              'byte_order="LittleEndian" header_type="UInt64"><UnstructuredGrid>'
              '<FieldData><Array type="String" Name="proj4" NumberOfTuples="1" format="ascii">'
              f'{projection}</Array></FieldData><Piece NumberOfPoints="{nv}" NumberOfCells="{n}">'
              f'<Points>{points}</Points><Cells>{cells}</Cells><CellData>{data}</CellData>'
              '</Piece></UnstructuredGrid><AppendedData encoding="raw">_')
    with open(path, 'wb') as out:
        out.write(header.encode())
        for dtype, size, blocks in arrays:
            out.write(struct.pack('<Q', size))
            written = 0
            for block in blocks():
                raw = np.asarray(block, dtype=dtype).tobytes(order='C')
                out.write(raw)
                written += len(raw)
            if written != size:
                raise ValueError('VTU array length changed during export')
        out.write(b'</AppendedData></VTKFile>')
