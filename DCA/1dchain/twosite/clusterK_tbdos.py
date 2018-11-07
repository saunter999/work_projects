#!/usr/bin/env python
from scipy import *
from pylab import *
import random as r

def K0_cluster(x):
	if x<0.25 or x>=0.75 :return True

def Kpi_cluster(x):
	if x>=0.25 and x<0.75 :return True


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
	K0dos=zeros(len(ebin));K0num=0.0
	Kpidos=zeros(len(ebin));Kpinum=0.0
	v=2.0*pi
	for i in range(M):
	    x=r.random()
	    k=x*v
	    kidx=0
	    if(K0_cluster(x)):   
		K0num+=1;kidx=1
	    if(Kpi_cluster(x)):  
		Kpinum+=1;kidx=2 
	    ek=-2.0*t*cos(k)
	    for ind,e in enumerate(ebin): 
	       if ek-e>=-De/2.0 and ek-e<De/2.0: 
		  if kidx==1:  K0dos[ind]+=1.0
		  if kidx==2:  Kpidos[ind]+=1.0 

	print "# of K0 cluster k points:",K0num
	print "# of Kpi cluster k points:",Kpinum
	print "their sum:",K0num+Kpinum

	plot( ebin,K0dos/(De*M),'k',lw=3,label="K=0")
	plot( ebin,Kpidos/(De*M),'g',lw=3,label="K=pi")
	plot( ebin,(K0dos+Kpidos)/(De*M),'m--',lw=1.5,label="sum" )
	xlim([-2.05*t,2.05*t])
	xlabel(r"$\omega$",size='x-large')
	ylabel("DOS",size='x-large')
	legend(loc=0.0)
	savefig("dos_clusterK.png")

def method2(M,Ne,emax):
	dlt=1.0e-2
	M=int(M);Ne=int(Ne)
	print "Lorentzian broadening method"
	print "# of Monte carlo steps is:",M
	print "# of engery bins sampled is:",Ne
	emesh=linspace(-emax,emax,Ne)
	v=2.0*pi
	AwK0=[]
	AwKpi=[]
	for omg in emesh:
	    GK0=0.0j
	    GKpi=0.0j
	    for i in range(M):
		x=r.random()
		k=x*v
		ek=-2.0*t*cos(k)
		if(K0_cluster(x)):   
		    GK0+=1./(omg-ek+1.0j*dlt)
		if(Kpi_cluster(x)):  
		    GKpi+=1./(omg-ek+1.0j*dlt)
	    GK0=GK0/M
	    GKpi=GKpi/M
	    AwK0.append(-1.0/pi*GK0.imag)
	    AwKpi.append(-1.0/pi*GKpi.imag)
	AwK0=array(AwK0)
	AwKpi=array(AwKpi)
	plot(emesh,AwK0,'k',lw=3,label="K=0")
	plot(emesh,AwKpi,'g',lw=3,label="K=pi")
	plot(emesh,AwK0+AwKpi,'m--',lw=1.5,label="sum")
	xlabel(r"$\omega$",size='x-large')
	ylabel("DOS",size='x-large')
	legend(loc=0.0)
	savefig("dos_clusterKMethod2.png")
	   
	




if __name__=="__main__":
	t=1.0;Ne=200;emax=3.*t
	M=1e5
#	clusterK_dos(M,Ne,emax)
	method2(M,Ne,emax)
	show()
