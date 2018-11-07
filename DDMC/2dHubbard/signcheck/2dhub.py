#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
from numpy import random
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


def G0rtau_gen(tauls):
	"""
	G0 expression in (r,tau) representation
	"""
	G0rtau=zeros((N*N,Ntau),dtype=complex)
	rls=[];ir=0
	for ix in arange(N):
	    for iy in arange(N):
		rls.append([ix,iy])
		for ik,ek in enumerate(ekls): 
		    G0rtau[ir,:]+=G0ktau(ek,tauls)*exp(1.0j*(ix*kls[ik][0]+iy*kls[ik][1]))  ##vector sum in python
		ir+=1
	G0rtau=G0rtau/(N-0.0)**2
	return G0rtau,rls
		    

def Grandpotential(n,ekls):
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
	   Ns=int(1e6)
	   Ne=int(8e6)
	   Nls=arange(Ns,Ne,Ne/6)	
	   print
	   print "calculating Grand potential correction at order "+str(n)
#	   omega=-U**2/2.0*SimpleMC(Nitt)
	   f=open("sign.txt",'w')
	   for Nitt in Nls:
	       print "Nitt=",Nitt
	       sign=MCMC(Nitt)
	       print>>f,Nitt,sign
	   
	return omega
		
def MCMC(Nitt):
	warm=5000
	measure=10
	sign=SampleSign(Nitt,warm,measure)
	return sign
	

def SampleSign(Nitt,warm,measure):
	sgn=0.0
	Naver=0.0
	itau=random.randint(0,Ntau)
	ix=random.randint(0,N);iy=random.randint(0,N)
#	tau=tauls[random.randint(0,Ntau)]
#	r=[random.randint(0,N),random.randint(0,N)]
	for itt in range(Nitt):
	    itaup=random.randint(0,Ntau)
	    ixp=random.randint(0,N);iyp=random.randint(0,N)
#	    taup=tauls[random.randint(0,Ntau)]
#	    rp=[random.randint(0,N),random.randint(0,N)]
	    P=min([abs(kernelMr_set(itaup,ixp*(N-1)+iyp))/abs(kernelMr_set(itau,ix*(N-1)+iy)),1])
	    if P>rand():
	       itau=itaup
	       ix=ixp;iy=iyp
	    if itt>warm and itt%measure==0:
		sgn+=np.sign(kernelMr_set(itau,ix*(N-1)+iy))
		Naver+=1.0
	return sgn/Naver
	
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

def kernelMr_set(itau,ir):
	Mr=zeros((4,4),dtype=complex)
	Mr[0,0]=-G0rtau[0,-1]
	Mr[1,1]=Mr[0,0]
	Mr[2,2]=Mr[0,0]
	Mr[3,3]=Mr[0,0]
	Mr[0,2]=G0rtau[ir,itau]
	Mr[1,3]=Mr[0,2]
	Mr[2,0]=-G0rtau[ir,Ntau-1-itau] 
	Mr[3,1]=Mr[2,0]
#	print r,tau,np.sign(LA.det(Mr)-n_nint**4)
	return LA.det(Mr).real-n_nint**4
	
	    
if __name__=="__main__":
	##model parameter
	N=20			   ##N--linear dimension of 2d square lattice
	t,U=1.0,2.0     	   ##t--hopping parameter;U--Hubbard on-site U.
#	beta,mu=100,0.0
	beta,mu=8,0.55978          ##beta--inverse temperature;mu--chemical potential.(mu=0.55978 gives density n=0.875 according to the paper.)
	Ntau=2000         	   ##Ntau--number of tau points in [0,beta)

	##Calculating G0 in (r,tau) space
	tauls=linspace(0,beta,Ntau)
	ekls,kls=ekdispersion(N)
	G0rtau,rls=G0rtau_gen(tauls) 
	n_nint=-G0rtau[0,-1].real
	
	print "Occupation number of H0 per spin is:"
	print n_nint

       	##Calculating Grand canonical potential using DDMC 
	Nexp=3			   ##Nexp--Largest expansion order
	omega=zeros(Nexp)          ##omega--stores the value of Grand canonical potential order by order
	for i in range(Nexp):
	    omega[i]=Grandpotential(i,ekls)
	print 
	print "Grand canonical potential order by order is "
	print omega
	print 
	print "Total canonical poetential is"
	print sum(omega)
	f=open("Omega_result.txt",'w')
	print>>f,omega
	print>>f, sum(omega)

	show()
	
	
