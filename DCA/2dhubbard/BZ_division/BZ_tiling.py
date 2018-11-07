#!/usr/bin/env python
from scipy import *
from pylab import *
import random as r

def K00_cluster(x1,x2):
	if x1<0.25 and (x2<0.25 or x2>0.75) :return True
	if x1>0.75 and (x2<0.25 or x2>0.75):return True
	
def Kpi0_cluster(x1,x2):
	if (x1>0.25 and x1<0.75) and (x2<0.25 or x2>0.75):return True

def K0pi_cluster(x1,x2):
	if (x1<0.25 or x1>0.75) and (x2>0.25 and x2<0.75):return True
	

def Kpipi_cluster(x1,x2):
	if (x1>0.25 and x1<0.75) and (x2>0.25 and x2<0.75):return True

if __name__=="__main__":
	M=10e4
	K00=[]
	Kpi0=[]
	K0pi=[]
	Kpipi=[]
	for i in range(int(M)):
	    x1=r.random()
	    x2=r.random()
	    if(K00_cluster(x1,x2)): 
	         K00.append([x1,x2])
	    if(Kpi0_cluster(x1,x2)): 
	         Kpi0.append([x1,x2])
	    if(K0pi_cluster(x1,x2)): 
	         K0pi.append([x1,x2])
	    if(Kpipi_cluster(x1,x2)): 
	         Kpipi.append([x1,x2])

	print "total number of k points:",int(M)
	print "sum of number of k points in four cluster momentum regions:"
	print len(K00)+len(Kpi0)+len(K0pi)+len(Kpipi)	
	K00=array(K00).transpose()
	Kpi0=array(Kpi0).transpose()
	K0pi=array(K0pi).transpose()
	Kpipi=array(Kpipi).transpose()
	plot(K00[0],K00[1],'ro')
	plot(Kpi0[0],Kpi0[1],'go')
	plot(K0pi[0],K0pi[1],'bo')
	plot(Kpipi[0],Kpipi[1],'mo')

	xlim([0,1])
	ylim([0,1])
	axvline(x=0.25,c='k',ls='--')
	axvline(x=0.75,c='k',ls='--')
	axhline(y=0.25,c='k',ls='--')
	axhline(y=0.75,c='k',ls='--')
	savefig("tilingBZ.png")
	show()
	
