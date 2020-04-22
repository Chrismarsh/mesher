#!/usr/bin/env python

import sys, os,pdb
import json
import numpy as np

# Set up CL arguments
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-o", "--outfile", required=False,
                help="File for output (default overwrites reference file).")

parser.add_argument("-ref", "--reffile", required=True,
                help="Reference file containing original mesh.")

parser.add_argument("-i", "--infile", required=True,
                help="File containing parameters to be added.")

import matplotlib as mpl
mpl.use('AGG')  # non-gui display (much faster)
import matplotlib.pyplot as plt


def append_new_param_to_file(args):
    """Read dictionary of arguments, find parameters, add new parameters to reference json file"""

    with open(args["infile"]) as f:
        param_in = json.load(f)
      
    with open(args["reffile"]) as f:
        param_ref = json.load(f)

    # Add parameters to param_ref
    param_ref.update(param_in)

    with open(args["outfile"],'w') as f:
        json.dump(param_ref, f, indent=4)


if __name__=="__main__":

    # Parse the input arguments
    args = vars(parser.parse_args())
    if (args["outfile"]==None):
        args["outfile"] = args["reffile"]

    append_new_param_to_file(args)
