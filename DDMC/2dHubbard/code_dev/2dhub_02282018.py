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
	G0rtau=zeros((N*N,Ntau),dtype=complex)
	rls=[];ir=0
	for ix in arange(N):
	    for iy in arange(N):
		rls.append([ix,iy])
		for it,tau in enumerate(tauls): 
		    for ik,ek in enumerate(ekls): 
		    	G0rtau[ir,it]+=G0ktau(ek,tau)*exp(1.0j*(ix*kls[ik][0]+iy*kls[ik][1]))
		ir+=1
#	G0rtau=zeros(Ntau,dtype=complex)
#	for iz,tau in enumerate(tauls): 
#	    for ik,ek in enumerate(ekls): 
#		G0rtau[iz]+=G0ktau(ek,tau)
	G0rtau=G0rtau/(N-0.0)**2
	return G0rtau,rls
		    
	
#def fermi(x):
#	cut=beta*8*t  
#	if x>cut : return 0.0
#	if x<-1*cut: return 1.0
#	return 1./(exp(x)+1.)

#def check(ekls):
#	n=0.0
#	for ek in ekls:
#	    n+=fermi(beta*(ek-mu))
#	return n/(N-0.0)**2
	    
if __name__=="__main__":
	##model parameter
	N=10  			   ##N--linear dimension of 2d square lattice
	t,U=1.0,2.0     	   ##t--hopping parameter;U--Hubbard on-site U.
	beta,mu=100,0.0
##	beta,mu=8,0.55978          ##beta--inverse temperature;mu--chemical potential.(mu=0.55978 gives density n=0.875 according to the paper.)
	Ntau=20         	   ##Ntau--number of tau points in [0,beta)

	tauls=linspace(0,beta,Ntau)
	ekls,kls=ekdispersion(N)
	G0rtau,rls=G0rtau_gen(tauls,ekls,kls)
#	print G0rtau[0,0,:].real
#	print G0rtau[0,0,:].imag
	plot(tauls,G0rtau[0,:].real,'o-')
	xlim([0,beta+10])
	axvline(x=beta,ls='--')
	print G0rtau[0,-1]
	show()
	
	
