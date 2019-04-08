#!/usr/bin/env python

import sys, os
import json
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.csgraph as graph

import matplotlib as mpl
mpl.use('AGG')  # non-gui display (much faster)
import matplotlib.pyplot as plt


def print_usage():
    print('USAGE: ')
    print('  permutation_tools.py <json_mesh_input_file> [json_mesh_output_file]')
    return

def append_local_cell_id_to_mesh_file(infile,outfile):
    """Read file, find a desired permutation of cell faces, add new permuted ids to json file"""

    with open(infile) as f:
        mesh = json.load(f)

    mesh['mesh']['cell_global_id'] = compute_minimum_bandwidth_permutation(mesh['mesh']['neigh'])

    with open(outfile,'w') as f:
        json.dump(mesh, f, indent=4)

    return


# permutation functions
def compute_minimum_bandwidth_permutation(neighbor_list):
    """Computes the permutation that minimizes the bandwidth of the connectivity matrix
    - uses reverse CutHill-McKee algorithm, which needs CSC or CSR sparse matrix format"""

    A = convert_neighbor_list_to_csr_matrix(neighbor_list)

    # Note we are always dealing with undirected (symmetric) graphs
    permutation = graph.reverse_cuthill_mckee(A,symmetric_mode=True)

    print_bandwidth_before_after_permutation(A,permutation)

    return permutation.tolist()

def convert_neighbor_list_to_csr_matrix(neighbor_list):
    """Get a scipy.sparse.csr_matrix from a list of neighbors"""

    N = len(neighbor_list)
    i = []
    j = []

    # easier to construct in coo format
    count = 0;
    for ii, neigh in enumerate(neighbor_list):
        for jj in neigh:
            if (jj !=-1):
                i.append(ii)
                j.append(jj)
                count=count+1
    val = np.ones(count)
    A = sparse.coo_matrix((val,(i,j)),shape=(N,N))

    return sparse.csr_matrix(A)


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

    if len(sys.argv) == 1:
        print_usage()
        exit(-1)

    infile = sys.argv[1]

    outfile = infile
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

    append_local_cell_id_to_mesh_file(infile, outfile)
