#!/usr/bin/env python

import sys, os
import json
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.csgraph as graph
import scipy.sparse.linalg as linalg
import importlib

# This is the path for the metis shared obj we detected during install
metis_config_path = os.path.dirname(os.path.abspath(__file__)) + '/metis-config.py'
metisconfig = importlib.machinery.SourceFileLoader('metisconfig', metis_config_path)
metisconfig = metisconfig.load_module()

if metisconfig.is_conan_build:
    filename = os.path.split(metisconfig.METIS_DLL)
    filename = filename[1] #the so name
    path = os.path.dirname(os.path.abspath(__file__)) + '/../lib/' + filename
    os.environ["METIS_DLL"] = path
else:
    os.environ["METIS_DLL"] = metisconfig.METIS_DLL
import metis

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

import matplotlib as mpl
mpl.use('AGG')  # non-gui display (much faster)
import matplotlib.pyplot as plt

# Version of the mesh that is produced
MESH_MAJOR = "1"
MESH_MINOR = "2"
MESH_PATCH = "0"


def append_global_cell_id_to_mesh_file(args):
    """Read dictionary of arguments, find a desired permutation of cell faces, add new permuted ids to json file"""

    with open(args["infile"]) as f:
        mesh = json.load(f)

    if args["type"]=="rcm":
        print(" Performing RCM bandwidth minimization:")
        mesh['mesh']['cell_global_id'] = compute_rcm_permutation(mesh['mesh']['neigh'])
    elif args["type"]=="metis":
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

    if nrank <= 1:
        raise Exception("Metis partitioning requires >1 partitions. Use -n / --nrank to set this")

    N = len(neighbor_list) # number of rows/columns

    trim_neighbor_list = [x[:] for x in neighbor_list]
    for nlst in trim_neighbor_list:
        while -1 in nlst:
            nlst.remove(-1)
    # partition by metis
    m = metis.part_graph(trim_neighbor_list,nrank)

    # add index to partitions sort, and extract new ordering
    L = [ (m[1][i],i) for i in range(len(m[1])) ]
    L.sort()
    sorted_l,permutation = zip(*L)

    # Get sizes of each partition
    rank_size = np.zeros(nrank,dtype=int)
    for rank in range(nrank):
        rank_size[rank] = m[1].count(rank)

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

############################################################################
############################################################################
if __name__=="__main__":

    # Parse the input arguments
    args = vars(parser.parse_args())
    if (args["outfile"]==None):
        args["outfile"] = args["infile"]

    append_global_cell_id_to_mesh_file(args)
