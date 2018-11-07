#!/usr/bin/env python
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()

if __name__ == '__main__':
    rank = comm.Get_rank()

    print size, rank
