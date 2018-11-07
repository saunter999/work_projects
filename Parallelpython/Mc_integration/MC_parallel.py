#!/usr/bin/env python
from scipy import *
from scipy import integrate
from mpi4py import MPI

def func(x):
        """function has to be postive definite in the range [0,beta].
        """
#        return exp(-x)
        return x 
#       return cosh(x)
#       return exp(x)

def SimpleMC(Nitt):
	if Parallel:
	    #print "Parallel implementation"
	    comm=MPI.COMM_WORLD
	    rank= comm.Get_rank()
	    size=comm.Get_size()
        sum=0.0

        for itt in range(Nitt):
             tau=rand()*beta
             ftau=func(tau)
             sum+=ftau

	###send all results to rank=0
	if Parallel:
	    #print "Gathering sum"
	    for i in range(1,size):
		comm.send(sum,dest=0)
	    if rank==0:
		for i in range(1,size):
		    sum+=comm.recv(source=i)
		print 'exact=',exact
	        print "Nitt=",Nitt*size
	        print sum/(Nitt*size)*beta

	if Parallel==False:
	        print "Nitt=",Nitt
	        print "MC=",sum/Nitt*beta


if __name__=="__main__":
	Parallel = True
	beta=8
	tls=linspace(0,beta,1000)
        yls=array([func(tau) for tau in tls])
        exact=integrate.simps(yls,tls)

	
	Nitt=1000
	SimpleMC(Nitt)
