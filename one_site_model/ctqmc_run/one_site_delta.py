#!/usr/bin/env python
from scipy import *
from pylab import *




if __name__=="__main__":
	t=1.0
	beta=50
	Nom=2000
	f=open("Delta.inp",'w')	
	for i in range(Nom):
	    om=(2*i+1)*pi/beta
	    delta=t**2/(1.0j*om) 
	    print>>f,om,delta.real,delta.imag
