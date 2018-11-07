#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import random





if __name__=="__main__":
	Nitt=10000
	Nup=0
	Ndn=0
	for it in range(Nitt):
	    if sign(rand()-0.5)==-1:
		Ndn+=1
	    else: Nup+=1

	N=6
	print random.randint(0,N)
	print Nup,Ndn,Nup+Ndn,Nitt
