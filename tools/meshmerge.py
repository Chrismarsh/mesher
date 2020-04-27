#!/usr/bin/env python
import json
from pathlib import Path

# Set up CL arguments
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-o", "--outfile", required=False,
                help="File for output. Default is an filename from all input files.")

parser.add_argument("infile",
                    nargs='+',
                    help="Files containing parameters to be merged.")


def append_new_param_to_file(args):
    """Read dictionary of arguments, find parameters, add new parameters to reference json file"""

    params = args["infile"]

    # load the first to append to
    with open(params[0]) as f:
        param_merged = json.load( f )

    for i in range(1, len(params)):
        with open(params[i]) as f:
            param_in = json.load(f)
            # Add parameters to param_ref
            param_merged.update(param_in)

    lengths = []

    print('Merged param file has:')
    for param, data in param_merged.items():

        l = len(data)
        lengths.append(l)
        print(F"{param}: {l}")

    if len(set(lengths)) > 1:
        print('ERROR: Not all parameters have the same length, which they must.')

    with open(args["outfile"],'w') as f:
        json.dump(param_merged, f, indent=4)


if __name__=="__main__":

    # Parse the input arguments
    args = vars(parser.parse_args())
    if args["outfile"] is None:
        name = ''
        for f in args["infile"]:
            filename = Path(f).resolve().stem
            name = name + filename + '+'
        args["outfile"] = '%s.param' % name[:-1]

    append_new_param_to_file(args)
