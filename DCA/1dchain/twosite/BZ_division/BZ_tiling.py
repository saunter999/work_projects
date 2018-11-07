#!/usr/bin/env python
from scipy import *
from pylab import *
import random as r

def K0_cluster(x):
	if x<0.25 or x>=0.75 :return True

def Kpi_cluster(x):
	if x>=0.25 and x<0.75 :return True

if __name__=="__main__":
	M=10e5
	K0=[]
	Kpi=[]
	for i in range(int(M)):
	    x=r.random()
	    if(K0_cluster(x)): 
	         K0.append([x])
	    if(Kpi_cluster(x)): 
	         Kpi.append([x])

	print "total number of k points:",int(M)
	NK0=len(K0)
	NKpi=len(Kpi)
	print "sum of number of k points in K=0 cluster momentum regions:"
	print NK0
	print "sum of number of k points in K=pi cluster momentum regions:"
	print NKpi	
	y=[0.001]*NK0
	plot(K0,y,'ro')
	y=[0.001]*NKpi
	plot(Kpi,y,'go')

	xlim([0,1])
#	ylim([0,1])
	axvline(x=0.25,c='k',ls='--')
	axvline(x=0.75,c='k',ls='--')
	savefig("tilingBZ.png")
	show()
	
