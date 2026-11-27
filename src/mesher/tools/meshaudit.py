"""Inspect a CHM mesh and explicitly rebuild missing/corrupt neighbours."""

import argparse
import json
import os
from pathlib import Path
import tempfile

from mesher.mesh_topology import audit_mesh, DEFAULT_CHUNK_SIZE


def _write_new_mesh(path, document):
    """Publish a complete new file, refusing to overwrite an existing output."""
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '-', dir=str(path.parent))
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(document, stream, indent=4, allow_nan=False)
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


# Reconstruct adjacency to audit a saved CHM mesh. A requested repair writes a separate file so
# the original remains available for diagnosis.
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', help='CHM .mesh JSON file')
    parser.add_argument('--rebuild-neighbours', action='store_true', help='Reconstruct neighbour rows without changing triangles or vertices')
    parser.add_argument('-o', '--outfile', help='New output file; existing files are never overwritten')
    parser.add_argument('--json', action='store_true', help='Print the audit report as JSON')
    parser.add_argument('--chunk-size', type=int, default=DEFAULT_CHUNK_SIZE, help='Maximum triangles per edge-sorting run')
    parser.add_argument('--boundary-tolerance', type=float, default=1e-8, help='Coincidence and hanging-node tolerance in mesh coordinate units (default: 1e-8)')
    args = parser.parse_args(argv)
    if bool(args.outfile) != args.rebuild_neighbours:
        parser.error('--rebuild-neighbours and --outfile must be used together')
    if args.chunk_size < 1:
        parser.error('--chunk-size must be positive')
    if args.outfile and Path(args.infile).resolve() == Path(args.outfile).resolve():
        parser.error('Use a separate output file; the input mesh is never overwritten')
    try:
        with open(args.infile) as stream:
            document = json.load(stream)
        report, neighbours = audit_mesh(document, args.chunk_size, args.boundary_tolerance)
        if args.rebuild_neighbours and neighbours is not None and not report['errors']:
            document['mesh']['neigh'] = neighbours.tolist()
            _write_new_mesh(args.outfile, document)
            report['previous_neighbour_errors'] = report['neighbour_errors']
            report['neighbour_errors'] = []
            report['valid'] = True
            report['output'] = args.outfile
            report.pop('recovery', None)
    except (OSError, TypeError, ValueError) as exc:
        parser.exit(1, f'meshaudit: {exc}\n')
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        for key in ('nvertex', 'nelem', 'boundary_edges', 'components', 'largest_component_sizes', 'hanging_nodes', 'coincident_boundary_vertices', 'output'):
            if key in report:
                print(f'{key}: {report[key]}')
        for message in report['warnings']:
            print(f'WARNING: {message}')
        for message in report['errors'] + report['neighbour_errors']:
            print(f'ERROR: {message}')
        if not report['valid'] and 'recovery' in report:
            print(report['recovery'])
        print('Domain coverage and arbitrary triangle overlaps are not checked.')
        print('Audit passed.' if report['valid'] else 'Audit failed; no repaired mesh was written.')
    return 0 if report['valid'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
