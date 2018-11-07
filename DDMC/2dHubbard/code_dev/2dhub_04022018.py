#!/usr/bin/env python
from scipy import *
from pylab import *
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


def G0rtau_gen(tauls,ekls,kls):
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
		    
def plotG0rtau(rs):
	plot(tauls,G0rtau[rs,:].real,'o-',label='real')
	plot(tauls,G0rtau[rs,:].imag,'*-',label='imag')
	xlim([0,beta+2])
	legend(loc=0)
	axvline(x=beta,ls='--')

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
	
	   
	return omega
		

	
	    
if __name__=="__main__":
	##model parameter
	N=10  			   ##N--linear dimension of 2d square lattice
	t,U=1.0,2.0     	   ##t--hopping parameter;U--Hubbard on-site U.
#	beta,mu=100,0.0
	beta,mu=8,0.55978          ##beta--inverse temperature;mu--chemical potential.(mu=0.55978 gives density n=0.875 according to the paper.)
	Ntau=20         	   ##Ntau--number of tau points in [0,beta)

	##Calculating G0 in (r,tau) space
	tauls=linspace(0,beta,Ntau)
	ekls,kls=ekdispersion(N)
	G0rtau,rls=G0rtau_gen(tauls,ekls,kls)
	n_nint=-G0rtau[0,-1].real
	print "Occupation number per spin is:"
	print n_nint

       	##Calculating Grand canonical potential using DDMC 
	Nexp=2			   ##Nexp--Largest expansion order
	omega=zeros(Nexp)          ##omega--stores the value of Grand canonical potential order by order
	for i in range(Nexp):
	    omega[i]=Grandpotential(i,ekls)
	print 
	print "Grand canonical potential order by order is "
	print omega
	print 
	print "Total canonical poetential is"
	print sum(omega)
#	plotG0rtau(0)
#	print G0rtau[0,-1]
	show()
	
	
