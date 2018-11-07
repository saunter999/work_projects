#!/usr/bin/env python
from mpi4py import MPI
comm=MPI.COMM_WORLD
rank= comm.Get_rank()
size=comm.Get_size()

#data=None
#if rank==0:
#	comm.send("helloworld",dest=1)
#if rank==1:
#	print 'on task',rank,'before recv:   data = ',data
#	data=comm.recv(source=0)
#	print 'on task',rank,'after recv:   data = ',data



data=rank
for i in range(1,size):
    comm.send(data,dest=0)

if rank==0:
	for i in range(1,size):
	    data+=comm.recv(source=i)
	print data
