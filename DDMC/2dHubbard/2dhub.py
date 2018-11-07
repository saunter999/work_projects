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

	if n>=2:
	   Ns=int(1e6)
	   Ne=int(8e7)
	   Nls=arange(Ns,Ne,Ne/4)	
	   print
	   print "calculating Grand potential correction at order "+str(n)
	   f=open("Mc_convergence.txt",'w')
	   print>>f,"n||Nitt||sign||norm||sign*norm||omega"
	   for Nitt in Nls:
	       print "Nitt=",Nitt
	       sign,norm=MCMC(n,Nitt)
	       print sign,norm,sign*norm,(-1)**(n+1)*U**n/n*sign*norm
	       print>>f,n,Nitt,sign,norm,sign*norm,(-1)**(n+1)*U**n/n*sign*norm
	   omega=(-1)**(n+1)*U**n/n*sign*norm
	   
	return omega
		
def MCMC(n,Nitt):
	warm=5000
	measure=20
	Dorder_init=0;c_norm=0.3
	sign=SampleSign(n,Nitt,warm,measure)
	norm=Samplenorm(n,Nitt,warm,measure,Dorder_init,c_norm)
	return sign,norm
	

def Samplenorm(n,Nitt,warm,measure,Dorder_init,c_norm):
	N0=0.0
        N1=0.0
        Dorder=Dorder_init
        Naver=0.0
        acc1=0.0
        acc2=0.0

        for itt in range(Nitt):
            if Dorder==0:
		   itau=random.randint(0,Ntau)
		   ix=random.randint(0,N);iy=random.randint(0,N)
                   P=min([beta*N**2*abs(kernelMr_set(n,itau,ix*(N-1)+iy))/c_norm,1])
                   if P>rand():#Metropolis
                    #accept:
                        Dorder=1
                        acc1+=1
                    #reject -> do nothing
            else:
                   P=min([c_norm/(beta*N**2*abs(kernelMr_set(n,itau,ix*(N-1)+iy))),1])
                   if P>rand():#Metropolis
                    #accept:
                        Dorder=0
                        acc2+=1
                    #reject -> do nothing
#           print acc1,acc2,itt

            if itt>warm and itt%measure==0:
                if Dorder==0:
                   N0+=1.
                if Dorder==1:
                   N1+=1.
                Naver+=1.0

	print "-----------------------------"
        print "number of N0||","number of N1||","sum of these two||"
        print N0,"(",N0/Naver,")",N1,"(",N1/Naver,")",Naver
        print "acceptance rate is:",acc1/Nitt,acc2/Nitt

	return c_norm*N1/N0

def SampleSign(n,Nitt,warm,measure):
	sgn=0.0
	Naver=0.0
	itau=random.randint(0,Ntau)
	ix=random.randint(0,N);iy=random.randint(0,N)
	for itt in range(Nitt):
	    itaup=random.randint(0,Ntau)
	    ixp=random.randint(0,N);iyp=random.randint(0,N)
	    P=min([abs(kernelMr_set(n,itaup,ixp*(N-1)+iyp))/abs(kernelMr_set(n,itau,ix*(N-1)+iy)),1])
	    if P>rand():
	       itau=itaup
	       ix=ixp;iy=iyp
	    if itt>warm and itt%measure==0:
		sgn+=np.sign(kernelMr_set(n,itau,ix*(N-1)+iy))
		Naver+=1.0
	return sgn/Naver


def kernelMr_set(n,itau,ir):
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
#	print abs(LA.det(Mr).real-n_nint**4)
	return LA.det(Mr).real-n_nint**4
	
	    
if __name__=="__main__":
	##model parameter
	N=25			   ##N--linear dimension of 2d square lattice
	t,U=1.0,2.0     	   ##t--hopping parameter;U--Hubbard on-site U.
#	beta,mu=100,0.0
	beta,mu=8,0.55978          ##beta--inverse temperature;mu--chemical potential.(mu=0.55978 gives density n=0.875 according to the paper.)
	Ntau=1200         	   ##Ntau--number of tau points in [0,beta)

	##Calculating G0 in (r,tau) space on discrete mesh
	print "Calculating G0 in (r,tau space) on discrete mesh"
	tauls=linspace(0,beta,Ntau)
	ekls,kls=ekdispersion(N)
	G0rtau,rls=G0rtau_gen(tauls) 
	n_nint=-G0rtau[0,-1].real
	print "Done~"
	
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
	
	
