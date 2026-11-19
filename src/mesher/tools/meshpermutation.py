#!/usr/bin/env python

import sys, os
import json
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.csgraph as graph
import scipy.sparse.linalg as linalg
import importlib.util

from mesher.mesh_topology import audit_mesh, MeshTopologyError, validate_neighbours

def _load_metis():
    # Validate inputs before loading METIS. Respect an explicitly selected DLL;
    # otherwise use the library found at build time, or the wrapper's discovery.
    config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'bin', 'metis-config.py')
    if 'METIS_DLL' not in os.environ and os.path.isfile(config_path):
        spec = importlib.util.spec_from_file_location('metisconfig', config_path)
        config = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(config)
        os.environ['METIS_DLL'] = config.METIS_DLL
    import metis
    return metis

# Set up CL arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", required=False,
                help="File for output (default overwrites input file).")

parser.add_argument("-t", "--type", required=False, choices={"rcm","nd","metis"},
                help="Type of reordering to perform.",
                nargs='?', const="rcm", type=str, default="rcm")

parser.add_argument("-i", "--infile", required=True,
                help="File for input.")

parser.add_argument("-n", "--nrank", required=False,
                help="Number of ranks to order for. (only works with metis order)",
                type=int, default=1)

# Version of the mesh that is produced
MESH_MAJOR = "1"
MESH_MINOR = "2"
MESH_PATCH = "0"


def append_global_cell_id_to_mesh_file(args):
    """Read dictionary of arguments, find a desired permutation of cell faces, add new permuted ids to json file"""
    global MESH_MAJOR, MESH_MINOR, MESH_PATCH

    with open(args["infile"]) as f:
        mesh = json.load(f)

    if args["type"]=="rcm":
        print(" Performing RCM bandwidth minimization:")
        mesh['mesh']['cell_global_id'] = compute_rcm_permutation(mesh['mesh']['neigh'])
    elif args["type"]=="metis":
        report, _ = audit_mesh(mesh, check_boundaries=False)
        errors = report['errors'] + report['neighbour_errors']
        if errors:
            raise MeshTopologyError('; '.join(errors) + '. Inspect the mesh with meshaudit before partitioning')
        for warning in report['warnings']:
            print('Warning: ' + warning)
        mesh['mesh']['local_size'], mesh['mesh']['cell_global_id'] = \
            compute_metis_permutation_N(mesh['mesh']['neigh'], args['nrank'])
        MESH_MAJOR = "2"
        MESH_MINOR = "0"
    elif args["type"]=="nd":
        print(" Performing ND fill-in minimization:")
        mesh['mesh']['cell_global_id'] = compute_nd_permutation(mesh['mesh']['neigh'])

    mesh['mesh']['partition_method'] = args["type"]
    mesh['mesh']['version'] = f'{MESH_MAJOR}.{MESH_MINOR}.{MESH_PATCH}'

    with open(args["outfile"],'w') as f:
        json.dump(mesh, f, indent=4)

    return

############################################################################
# permutation functions
def compute_rcm_permutation(neighbor_list):
    """Computes the permutation that minimizes the bandwidth of the connectivity matrix
    - uses Reverse CutHill-McKee (RCM) algorithm, which needs CSC or CSR sparse matrix format"""

    A = convert_neighbor_list_to_compressed_matrix(neighbor_list, sparse.csr_matrix)

    # Note we are always dealing with undirected (symmetric) graphs
    permutation = graph.reverse_cuthill_mckee(A,symmetric_mode=True)

    print_bandwidth_before_after_permutation(A,permutation)

    return permutation.tolist()

def compute_nd_permutation(neighbor_list):
    """Computes a nested dissection permutation (minimizes fill-in of matrix factors)"""
    A = convert_neighbor_list_to_compressed_matrix(neighbor_list,sparse.csc_matrix)

    # permutation is determined during the factorization
    LUperm = linalg.splu(A)

    # Not really necessary, but curious
    print_bandwidth_before_after_permutation(
        convert_neighbor_list_to_compressed_matrix(neighbor_list, sparse.csr_matrix),
        LUperm.perm_c )

    print("  Non-zeros in factor L: " + str(LUperm.L.nnz))

    return LUperm.perm_c.tolist()

def compute_metis_permutation_N(neighbor_list,nrank):
    """Computes an nrank partition using METIS
    - Uses metis deffault permutation for partitions"""

    neighbours = validate_neighbours(neighbor_list)
    N = len(neighbours)
    if not isinstance(nrank, (int, np.integer)) or isinstance(nrank, bool) or not 2 <= nrank <= N:
        raise MeshTopologyError(f'METIS requires 2 <= nrank <= number of triangles ({N}); got {nrank}')

    if nrank == N:
        # There is exactly one possible partition size: one face per rank.
        # Some METIS versions return empty parts for this tiny-graph case.
        labels = np.arange(N)
    else:
        trim_neighbor_list = [row[row >= 0].tolist() for row in neighbours]
        metis = _load_metis()
        m = metis.part_graph(trim_neighbor_list,nrank)
        labels = np.asarray(m[1])
    if labels.shape != (N,) or labels.dtype.kind not in 'iu':
        raise MeshTopologyError(f'METIS returned an invalid partition array; expected {N} integer labels')
    if labels.min() < 0 or labels.max() >= nrank:
        raise MeshTopologyError('METIS returned partition labels outside the requested rank range')
    rank_size = np.bincount(labels.astype(np.int64), minlength=nrank)
    if np.any(rank_size == 0):
        raise MeshTopologyError('METIS returned an empty partition; reduce nrank or inspect disconnected components')
    # This is new-position -> original-face-index, as consumed by CHM.
    permutation = np.argsort(labels, kind='stable').tolist()

    A = convert_neighbor_list_to_compressed_matrix(neighbor_list, sparse.csr_matrix)
    print_bandwidth_before_after_permutation(A,permutation)

    return rank_size.tolist(), permutation


def convert_neighbor_list_to_compressed_matrix(neighbor_list,compressed_format):
    """Get a sparse compressed_format from a list of neighbors (assuming aij input)
    - compressed_format is one of sparse.csr_matrix or sparse.csc_matrix"""

    N = len(neighbor_list)
    i = []
    j = []
    val = []

    # easier to construct in coo format
    count = 0;
    for ii, neigh in enumerate(neighbor_list):
        for jj in neigh:
            if (jj !=-1):
                i.append(ii)
                j.append(jj)
                count=count+1
                val.append( (count % 2) + 1 ) # entries filled with 1s and 2s
    A = sparse.coo_matrix((val,(i,j)),shape=(N,N))

    return compressed_format(A)


## Bandwidth comparison
def print_bandwidth_before_after_permutation(A,permutation):
    """Display bandwidth before and after permutation"""
    orig_band = compute_bandwidth_CSR(A)
    print("  Original bandwidth: ", orig_band)
    plot_mat_connectivity(A,"orig_order.png")
    Aperm = A[permutation,:][:,permutation]

    post_band = compute_bandwidth_CSR(Aperm)
    plot_mat_connectivity(Aperm,"post_order.png")
    print("  Permuted bandwidth: ", post_band)
    return

def compute_bandwidth_CSR(A):
    """Compute the maximum bandwidth from a CSR sparse matrix A"""
    assert(sparse.isspmatrix_csr(A))
    indices = A.indices
    indptr = A.indptr
    bandwidth = 0
    for row in range(A.shape[0]):
        row_indices = indices[indptr[row]:indptr[row+1]]
        if not len(row_indices):
            continue
        maxind = np.max(row_indices)
        if maxind-row > bandwidth:
            bandwidth = maxind-row
    return bandwidth

def plot_mat_connectivity(A,filename,**kwargs):
    """Plot the spy view of the connectivity matrix
    - TODO kwargs can be used to get values to pass into matplotlib functions"""

    import matplotlib as mpl
    mpl.use('AGG')
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(18,16), dpi= 80)
    plt.spy(A,markersize=0.05)
    Nnz = A.nnz
    Nrow = A.shape[0]
    ratio = Nnz/(Nrow*Nrow)
    plt.title("Nearest neighbour connectivity")
    plt.xlabel("Number of non-zeros: " + str(Nnz) + " (%.3f %%)" %(ratio*100))
    plt.savefig(filename)
    plt.close(fig)

############################################################################
############################################################################


def main():
    # Parse the input arguments
    args = vars(parser.parse_args())
    if args["outfile"] is None:
        args["outfile"] = args["infile"]

    try:
        append_global_cell_id_to_mesh_file(args)
    except MeshTopologyError as exc:
        parser.exit(2, f'meshpermutation: {exc}\n')


if __name__ == "__main__":
    main()
