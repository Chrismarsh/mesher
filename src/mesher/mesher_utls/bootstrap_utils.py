import os
import sys


def ensure_mesher_on_path():
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))


def check_worker_count(actual, expected, setting='MPI_nworkers'):
    if expected is not None and actual != int(expected):
        raise ValueError(f'MPI launcher started {actual} ranks, but {setting}={expected}')
