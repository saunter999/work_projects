#!/usr/bin/env python
from scipy import *
from numpy import linalg as LA
import numpy as np
from pylab import *




def Ham_diag(mu):
	H=zeros((8,8))
	H[1,1]=-mu;H[1,3]=t
	H[2,2]=-mu;H[2,3]=t
	H[3,1]=t;H[3,2]=t;H[3,3]=-mu*0
	H[4,4]=U-2.*mu;H[4,5]=t;H[4,6]=-t
	H[5,4]=t;H[5,5]=-mu*1
	H[6,4]=-t;H[6,6]=-mu*1
	H[7,7]=U-2.*mu
	Eig,v=LA.eig(H)	
        #idx=Eig.argsort()
        #Eig=Eig[idx]
        #v=v[:,idx]
        return (Eig,v.transpose())

def UnitaryM(v):
	Um=zeros((8,8))  
	for i in range(8):
	    Um[i,:]=v[i,:]
	return Um
	
def d_matrix():
	dm=zeros((8,8))
	dm[1,0]=1.0
	dm[4,2]=1.0
	dm[5,3]=1.0
	dm[7,6]=1.0
	return dm

def t_matrix(N):
	tm=zeros((N,N))
	for k in range(N):
	    for i in range(N):
	        for j in range(N):
		    for m in range(N):
		       tm[i,k]+=Um[i,j]*dm[j,m]*Um.transpose()[m,k]
	return tm
if __name__=="__main__":
	beta=50
	t=1.0;Uls=array([0.0,8.0]);
	muls=Uls/2.0
	omls=[];G0om=[];GUom=[];Nom=2000
	dm=d_matrix()

	for i in range(Nom):
	    omls.append( (2*i+1)*pi/beta )
	
	for idx,U in enumerate(Uls):
	     print U
	     mu=muls[idx]
	     Eig,v=Ham_diag(mu)
	     Um=UnitaryM(v)
	     expE=[];Z=0.0
	     for e in Eig:
		expE.append(exp(-beta*e))
	     expE=array(expE)
	     Z=sum(expE)
	     N=len(Eig)
	     tm=t_matrix(N)
 		
	     Gom=[]
	     for om in omls: 
		G=0.0j
		for i in range(N):
		    for j in range(N):
			if tm[i,j]!=0:
			     G+=tm[i,j]**2*(expE[i]+expE[j])/(1.0j*om+Eig[j]-Eig[i])
		Gom.append(G/Z)
	     Gom=array(Gom)
	     if U==0:
	       G0om=Gom
	       G0om=array(G0om)
	     else:
	       GUom=Gom
	       GUom=array(GUom)

	figure(1)
#	plot(omls,G0om.real,'o-',lw=2,label="real"+"U="+str(0))
#	plot(omls,G0om.imag,'^-',lw=2,label="imag"+"U="+str(0))
#	plot(omls,GUom.real,'o-',lw=2,label="real"+"U="+str(Uls[1]))
	plot(omls,GUom.imag,'^-',lw=2,label="imag"+"U="+str(Uls[1]))
	xlim([0,omls[-1]/20.])
	xlabel(r'$i\omega$',size='x-large')
	ylabel(r'$G(i\omega)$',size='x-large')
	legend(loc=0.0)
	fout=open("Gf.out",'w')
	for i,om in enumerate(omls):
	    print>>fout,om,G0om[i].real,G0om[i].imag
	fout.close()


	Sigma=[]
	fout=open("Sig.out",'w')
	for i,om in enumerate(omls):
	     sig=1.0/G0om[i]-1.0/GUom[i]
	     print>>fout,om,sig.real,sig.imag 
	     Sigma.append(sig)
	fout.close()
	Sigma=array(Sigma)
	figure(2)
	ref=loadtxt("/work/qhan/projects/one_site_model/ctqmc_run/ref_selfenergy.dat").transpose()
#	title(r'$\mu=$'+str(mu)+"U="+str(Uls[1]),size='x-large')
	plot(omls,Sigma.real,'o-',lw=2,label="real")
	plot(omls,Sigma.imag,'^-',lw=2,label="imag")
	plot(ref[0],ref[1],'o',label="ref_imag")
	#xlim([0,omls[-1]/20.])
	xlim([0,60])
	xlabel(r'$i\omega$',size='x-large')
	ylabel(r'$\Sigma(i\omega)$',size='x-large')
	legend(loc=0.0)
	show()
	
