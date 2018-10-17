#!/usr/bin/env python

import sys, os
import json
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.csgraph as graph

def print_usage():
    print('USAGE: ')
    print('  permutation_tools.py <json_mesh_input_file> [json_mesh_output_file]')
    return

def append_local_cell_id_to_mesh_file(infile,outfile):
    """Read file, find a desired permutation of cell faces, add new permuted ids to json file"""

    with open(infile) as f:
        mesh = json.load(f)

    mesh['mesh']['cell_local_id'] = compute_minimum_bandwidth_permutation(mesh['mesh']['neigh'])

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


if __name__=="__main__":

    if len(sys.argv) == 1:
        print_usage()
        exit(-1)

    infile = sys.argv[1]

    outfile = infile
    if len(sys.argv) == 3:
        outfile = sys.argv[2]

    append_local_cell_id_to_mesh_file(infile, outfile)
