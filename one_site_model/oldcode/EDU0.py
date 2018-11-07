#!/usr/bin/env python
from scipy import *
from numpy import linalg as LA
import numpy as np
from pylab import *



def G_anaexp(om):
	return 0.5/(1.0j*om+mu)+0.25/(1.0j*om+mu-sqrt(2)*t)+0.25/(1.0j*om+mu+sqrt(2)*t)

def Ham_diag():
	H=zeros((8,8))
	H[1,1]=-mu;H[1,3]=t
	H[2,2]=-mu;H[2,3]=t
	H[3,1]=t;H[3,2]=t;H[3,3]=-mu
	H[4,4]=U-2.*mu;H[4,5]=t;H[4,6]=-t
	H[5,4]=t;H[5,5]=-2.*mu
	H[6,4]=-t;H[6,6]=-2.*mu
	H[7,7]=U-3.*mu
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
	t=1.0;U=3;
	mu=U/2.0
	omls=[];Gom=[];Nom=2000
	for i in range(Nom):
	    omls.append( (2*i+1)*pi/beta )
	Eig,v=Ham_diag()
	Um=UnitaryM(v)
	dm=d_matrix()

#	for i in range(8):
#	    print Eig[i]
#	    print v[i,:]
#	print Um
	expE=[]
	for e in Eig:
	    expE.append(exp(-beta*e))
	expE=array(expE)
	Z=sum(expE)


	Gomexp=[]
	N=len(Eig)
	tm=t_matrix(N)
 		
	for om in omls: 
	  G=0.0j
	  for i in range(N):
	      for j in range(N):
		  if tm[i,j]!=0:
		       G+=tm[i,j]**2*(expE[i]+expE[j])/(1.0j*om+Eig[j]-Eig[i])
	  Gom.append(G/Z)
	  Gomexp.append(G_anaexp(om))
	Gom=array(Gom)
	Gomexp=array(Gomexp)
	title(r'$\mu=$'+str(mu),size='x-large')
	plot(omls,Gom.real,lw=2,label="real")
	plot(omls,Gomexp.real,'--',lw=3,label="real_exact")
	plot(omls,Gom.imag,lw=2,label="imag")
	plot(omls,Gomexp.imag,'--',lw=3,label="imag_exact")
	xlim([0,omls[-1]/20.])
	xlabel(r'$i\omega$',size='x-large')
	ylabel(r'$G(i\omega)$',size='x-large')
	legend(loc=0.0)
	fout=open("Gf.out",'w')
	for i,om in enumerate(omls):
	    print>>fout,om,Gom[i].real,Gom[i].imag
	fout.close()
	show()
	
