import os
import sys


# A launcher can execute a worker by filename from another directory. Add the package root so it
# imports the same Mesher checkout.
def ensure_mesher_on_path():
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))


# Fail before splitting work if the launcher supplied a different number of ranks from the one
# the parent used to plan the stage.
def check_worker_count(actual, expected, setting='MPI_nworkers'):
    if expected is not None and actual != int(expected):
        raise ValueError(f'MPI launcher started {actual} ranks, but {setting}={expected}')
