#!/usr/bin/env python
from scipy import *
from pylab import *





if __name__=="__main__":
	Nitt=10000
	Nup=0
	Ndn=0
	for it in range(Nitt):
	    if sign(rand()-0.5)==-1:
		Ndn+=1
	    else: Nup+=1

	print Nup,Ndn,Nup+Ndn,Nitt
