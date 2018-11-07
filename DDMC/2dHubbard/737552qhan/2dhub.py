#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
import numpy as np
"""
Benchmark calculation for PRL 119,045701(2017) Riccardo Rossi
"""
def ekdispersion(N):
	kmesh=linspace(0,2*pi,N,endpoint=False)
	ek=[]
	kls=[]
	for kx in kmesh:
	    for ky in kmesh:
		ek.append(-2.0*t*(cos(kx)+cos(ky)) )
		kls.append([kx,ky])
	ek=array(ek)
	return ek,kls


def G0ktau(ek,tau):
	"""
	G0 expression in (k,tau) representation
	"""
	x=ek-mu
	if x>0 : return -exp(-x*tau)/(1.0+exp(-x*beta))
	else   : return -exp(-x*(tau-beta))/(1.0+exp(x*beta))


def G0rtau_gen(tau,r):
	"""
	G0 expression in (r,tau) representation
	"""
	G0rtau=0.0
	for ik,ek in enumerate(ekls): 
	    G0rtau+=G0ktau(ek,tau)*exp(1.0j*(r[0]*kls[ik][0]+r[1]*kls[ik][1]))  ##vector sum in python
	G0rtau=G0rtau/(N-0.0)**2
	return G0rtau
		    

def Grandpotential(n,ekls,Nitt):
	omega=0.0
	if n==0:
	   print
	   print "calculating Grand potential for H0..."
	   for ek in ekls:
	       x=ek-mu
	       if x>0:
		  omega+=log(1+exp(-beta*x))
	       else :
		  omega+=-beta*x+log(1+exp(beta*x))
	   omega*=-2.0/beta
	   omega*=1./(N-0.0)**2
	if n==1:
	   print
	   print "calculating Grand potential correction at order "+str(n)
	   omega=U*n_nint*n_nint

	if n==2:
	   print
	   print "calculating Grand potential correction at order "+str(n)
	   omega=-U**2/2.0*SimpleMC(Nitt)
	   
	return omega
		
def SimpleMC(Nitt):
	sum=0.0
        for itt in range(Nitt):
             tau=rand()*beta
	     Mtau=0.0
#	     print itt,tau,kernelMr_set(tau,[1,0]),G0rtau_gen(tau,[1,0])
	     for ix in range(N):
	     	for iy in range(N):
		   Mtau+=kernelMr_set(tau,[ix,iy])
	     print itt,tau,Mtau.real
             sum+=Mtau.real
	print sum/Nitt*beta
        return sum/Nitt*beta

def kernelMr_set(tau,r):
	Mr=zeros((4,4),dtype=complex)
	Mr[0,0]=-G0rtau_gen(beta,[0,0])
	Mr[1,1]=Mr[0,0]
	Mr[2,2]=Mr[0,0]
	Mr[3,3]=Mr[0,0]
	Mr[0,2]=G0rtau_gen(tau,r)
	Mr[1,3]=G0rtau_gen(tau,r)
	Mr[2,0]=-G0rtau_gen(beta-tau,r)  
	Mr[3,1]=-G0rtau_gen(beta-tau,r)  
#	print r,tau,np.sign(LA.det(Mr)-n_nint**4)
	return LA.det(Mr)-n_nint**4
	
	    
if __name__=="__main__":
	##model parameter
	N=6 			   ##N--linear dimension of 2d square lattice
	t,U=1.0,2.0     	   ##t--hopping parameter;U--Hubbard on-site U.
#	beta,mu=100,0.0
	beta,mu=8,0.55978          ##beta--inverse temperature;mu--chemical potential.(mu=0.55978 gives density n=0.875 according to the paper.)
#	Ntau=20         	   ##Ntau--number of tau points in [0,beta)

	##Calculating G0 in (r,tau) space
#	tauls=linspace(0,beta,Ntau)
	ekls,kls=ekdispersion(N)
	G0r0=G0rtau_gen(beta,[0,0]) 
	n_nint=-G0r0.real
	
	print "Occupation number per spin is:"
	print n_nint

       	##Calculating Grand canonical potential using DDMC 
	Nitt=int(1e3)
	Nexp=3			   ##Nexp--Largest expansion order
	omega=zeros(Nexp)          ##omega--stores the value of Grand canonical potential order by order
	for i in range(Nexp):
	    omega[i]=Grandpotential(i,ekls,Nitt)
	print 
	print "Grand canonical potential order by order is "
	print omega
	print 
	print "Total canonical poetential is"
	print sum(omega)
	f=open("Omega_result.txt",'w')
	print>>f,omega
	print>>f, sum(omega)
#	plotG0rtau(0)
#	print G0rtau[0,-1]
	show()
	
	
