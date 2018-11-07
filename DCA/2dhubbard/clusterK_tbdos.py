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

def clusterK_dos(M,Ne,emax):
	M=int(M);Ne=int(Ne)
	print "Monte carlo sampling method"
	print "# of Monte carlo steps is:",M
	print "# of engery bins sampled is:",Ne
	emesh=linspace(-emax,emax,Ne)
	De=emesh[1]-emesh[0]
	ebin=[] ## this stores the midpoint of each energy bin.
	for ind,e in enumerate(emesh[:-1]):
	    ebin.append( (e+emesh[ind+1])/2.0 )
	ebin=array(ebin)
	K00dos=zeros(len(ebin));K00num=0.0
	Kpi0dos=zeros(len(ebin));Kpi0num=0.0
	K0pidos=zeros(len(ebin));K0pinum=0.0
	Kpipidos=zeros(len(ebin));Kpipinum=0.0
	v1=array([2.0*pi,0])	
	v2=array([0.,2.0*pi])	
	for i in range(M):
	    x1=r.random()
	    x2=r.random()
	    k=x1*v1+x2*v2
	    kidx=0
	    if(K00_cluster(x1,x2)):   
		K00num+=1;kidx=1
	    if(Kpi0_cluster(x1,x2)):  
		Kpi0num+=1;kidx=2 
	    if(K0pi_cluster(x1,x2)):  
		K0pinum+=1;kidx=3
	    if(Kpipi_cluster(x1,x2)): 
		Kpipinum+=1;kidx=4 
	    ek=-2.0*t*(cos(k[0])+cos(k[1]))
	    for ind,e in enumerate(ebin): 
	       if ek-e>=-De/2.0 and ek-e<De/2.0: 
		  if kidx==1:  K00dos[ind]+=1.0
		  if kidx==2:  Kpi0dos[ind]+=1.0 
		  if kidx==3:  K0pidos[ind]+=1.0
		  if kidx==4:  Kpipidos[ind]+=1.0  

	print "# of K00 cluster k points:",K00num
	print "# of Kpi0 cluster k points:",Kpi0num
	print "# of K0pi cluster k points:",K0pinum
	print "# of Kpipi cluster k points:",Kpipinum
	print "their sum:",K00num+Kpi0num+K0pinum+Kpipinum

	plot( ebin,K00dos/(De*M),'k',label="K=(0,0)")
	plot( ebin,Kpi0dos/(De*M),'r',label="K=(pi,0)")
	plot( ebin,K0pidos/(De*M),'b--',label="K=(0,pi)")
	plot( ebin,Kpipidos/(De*M),'g',label="K=(pi,pi)" )
	plot( ebin,(K00dos+Kpi0dos+K0pidos+Kpipidos)/(De*M),c='m',label="sum" )
	xlim([-4.01*t,4.01*t])
	xlabel(r"$\omega$",size='x-large')
	ylabel("DOS",size='x-large')
	legend(loc=0.0)
	savefig("dos_clusterK.png")




if __name__=="__main__":
	t=1.0;Ne=100;emax=8.*t
	M=3e6
	clusterK_dos(M,Ne,emax)
	show()
