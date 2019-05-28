#!/usr/bin/env python

import sys, os
import json
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.csgraph as graph
import scipy.sparse.linalg as linalg

# Set up CL arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", required=False,
                help="File for output (default overwrites input file).")

parser.add_argument("-t", "--type", required=False, choices={"rcm","nd"},
                help="Type of reordering to perform.",
                nargs='?', const="rcm", type=str, default="rcm")

parser.add_argument("-i", "--infile", required=True,
                help="File for input.")

import matplotlib as mpl
mpl.use('AGG')  # non-gui display (much faster)
import matplotlib.pyplot as plt

def append_global_cell_id_to_mesh_file(args):
    """Read dictionary of arguments, find a desired permutation of cell faces, add new permuted ids to json file"""

    with open(args["infile"]) as f:
        mesh = json.load(f)

    if args["type"]=="rcm":
        print(" Performing RCM bandwidth minimization:")
        mesh['mesh']['cell_global_id'] = compute_rcm_permutation(mesh['mesh']['neigh'])
    elif args["type"]=="nd":
        print(" Performing ND fill-in minimization:")
        mesh['mesh']['cell_global_id'] = compute_nd_permutation(mesh['mesh']['neigh'])

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

############################################################################

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
        maxind = np.max(indices[indptr[row]:indptr[row+1]])
        if maxind-row > bandwidth:
            bandwidth = maxind-row
    return bandwidth

def plot_mat_connectivity(A,filename,**kwargs):
    """Plot the spy view of the connectivity matrix
    - TODO kwargs can be used to get values to pass into matplotlib functions"""

    fig=plt.figure(figsize=(18,16), dpi= 80)
    plt.spy(A,markersize=0.05)
    Nnz = A.nnz
    Nrow = A.shape[0]
    ratio = Nnz/(Nrow*Nrow)
    plt.title("Nearest neighbour connectivity")
    plt.xlabel("Number of non-zeros: " + str(Nnz) + " (%.3f %%)" %(ratio*100))
    plt.savefig(filename)

if __name__=="__main__":

    # Parse the input arguments
    args = vars(parser.parse_args())
    if (args["outfile"]==None):
        args["outfile"] = args["infile"]

    append_global_cell_id_to_mesh_file(args)
