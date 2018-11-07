#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
import numpy as np

def safe_exp(x):
	if x<-50.: return 0.0
	return exp(x)

def fermionkernM(lamd,x,y):
	xdim=len(x);ydim=len(y)
	kernM=zeros((xdim,ydim))
	for i in range(xdim):
	    for j in range(ydim):
		kernM[i,j]=0.5*safe_exp(-lamd/2.0*x[i]*y[j])/cosh(lamd/2.0*y[j])
	return kernM

def svd_de(M):
	U,s,V=LA.svd(M)
	print "shpae of U,s,V:"
	print U.shape,s.shape,V.shape
	return U,s,V

def ortho_check(v1,v2):
	print "orthonormality check:"
	print np.dot(v1,v2)

def Gl_eva(U,Gtau,Dtau):
	Lmax=50
	Ntau=len(Gtau)
	Gl=zeros(Lmax)
	for i in range(Lmax):
	    for j in range(Ntau):
		Gl[i]+=U[j,i]*Gtau[j]
	return Gl
		
	     
def Gtau_fromGl(U,Gl,tauls):
	Lcut=19
	print "svd_lmax=",Lcut,"used to reconstruct G(tau)" 
	Gtau_gl=[]
	Ntau=len(tauls)
	for i in range(Ntau):
	    sumGtau=0.0
	    for j in range(Lcut):
	        sumGtau+=Gl[j]*U[i,j]
	    Gtau_gl.append(sumGtau)
	Gtau_gl=array(Gtau_gl)
	return Gtau_gl
	
	
if __name__=="__main__":
	Lamd=500
#	xmesh=linspace(-1,1,Nx)
	ftau=loadtxt("Gtau_qh.dat").transpose()
	beta=ftau[0][-1]
	Dtau=ftau[0][1]-ftau[0][0]
	tauls=ftau[0]
	Gtau=ftau[1]
	xmesh=[]
	for tau in ftau[0]:
	    xmesh.append(2.0*tau/beta-1.)
	xmesh=array(xmesh)
	Dx=xmesh[1]-xmesh[0]
	print "# of points of xmesh is",len(xmesh)
	Ny=len(xmesh)/2
	print "# of points of ymesh is",Ny
	ymesh=linspace(-1,1,Ny)
	fkernM=fermionkernM(Lamd,xmesh,ymesh)
	figure(0)
	imshow(fkernM.transpose(), interpolation='nearest',origin='lower') ##imshow somehow set x dim in y direction, so we plot matrix.transpose() instead.

	U,s,V=svd_de(fkernM)
	ortho_check(U[0,:],U[3,:])
	fs=open("singular_value.dat",'w')
	for i in range(len(s)):
	    print>>fs,s[i]
	fs.close()

	figure(1)
	semilogy(range(len(s)),s/s[0],label='singular_value')
	xlabel('L')
	ylabel('s_L/s_(L=0)')
	legend(loc=0)
	Lmax=6

	figure(2)
	for i in range(Lmax): 
  	   plot(xmesh,U[:,i]/U[-1,i],label="l="+str(i))
	xlabel('x')
	ylabel('u_L(x)/u_L(1)')
	legend(loc=0)
	figure(3)
	for i in range(Lmax): 
	   plot(ymesh,V[i,:]/V[i,-1],label="l="+str(i))
	xlabel('y')
	ylabel('v_L(y)/v_L(1)')
	legend(loc=0)
	
	Gl=Gl_eva(U,Gtau,Dtau)
	figure(4)
	plot(range(len(Gl)),Gl/Gl[0])
	xlabel('l')
	ylabel('Gl')
	print "Gl is:"
	print Gl
	Gtau_gl=Gtau_fromGl(U,Gl,tauls)
	figure(5)
	semilogy(tauls,-Gtau,'-',lw=2,label="input")
	semilogy(tauls,-Gtau_gl,'--',lw=3,label="IR")
	legend(loc=0)
	xlim([0,beta/10])
	xlabel(r'$\tau$')
	ylabel(r'$-G(\tau)$')
	
	show()
	
	
