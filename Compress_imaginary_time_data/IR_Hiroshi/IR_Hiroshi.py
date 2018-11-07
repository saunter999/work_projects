#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
import numpy as np

class Gtau:
	"""
	Class Gtau is to package G(tau) in a single class object  
	Description of its attribute: 
	NO inheriated attribute
	c.data(): 2d array containing of tau and G(tau)
	c.beta(): inverse temperature 
	c.tau() : 1d array containing imaginary time points tau  
	c.gt()  : 1d array containing Green's function G(tau) 
	c.Ntau(): length of number of tau points
	"""
	def __init__(self,Gdata):	
	    self.data=Gdata

	def beta(self):
	    return self.data[0][-1]

	def tau(self):
	    return self.data[0]

	def gt(self):
	    return self.data[1]

	def Ntau(self):
	    return self.data.shape[1]

		
class fkernel_linearmesh(Gtau):
	"""
	Class fkernel_linearmesh is to generate kernel K(x,y) on linearmesh of x and y. 
	Description of its attribute: 
	Inheritated attribute from Gtau
	c.fkernel(): return value--K(x,y)
	"""
	def __init__(self,Gdata):
	    Gtau.__init__(self,Gdata)

	def fkernel(self,lmbda,dimxyratio): 
	    beta=self.beta()
	    tau=self.tau()
	    Ntau=self.Ntau()
	    xmesh=[]
	    for t in tau:
		xmesh.append(2.0*t/beta-1.0)
	    xmesh=array(xmesh)
	    Ny=int(Ntau/dimxyratio)
	    ymesh=linspace(-1.,1.,Ny)
	    kernM=zeros((Ntau,Ny))
	    for i,x in enumerate(xmesh):
		for j,y in enumerate(ymesh):
		    kernM[i,j]=0.5*safe_exp(-lmbda/2.0*x*y)/cosh(lmbda/2.0*y)
#		    kernM[i,j]=0.5*exp(-lmbda/2.0*x*y)/cosh(lmbda/2.0*y)
	    print "kernel's shape is ", kernM.shape
	    figure(0)
	    imshow(kernM.transpose(), interpolation='nearest',origin='lower') ##imshow somehow set x dim in y direction, so we plot matrix.transpose() instead.
	    return kernM
		
class svd_de:
	"""
	Class svd_de is to do singular value decomposition of K(x,y)
	Description of its attribute: 
	c.svd_eva():return value--Matrix of U,s,V in the equation K=UsV^+.
	"""
	def __init__(self,fk):
	   self.fkm=fk
	   self.s=None
	   self.U=None
	   self.V=None

	def svd_eva(self,Lmax):
	    U,s,V=LA.svd(self.fkm)
	    print "shpae of U,s,V:"
	    print U.shape,s.shape,V.shape
	    self.s=s
	    self.U=U
	    self.V=V
	    figure(1)
	    xs=range(len(s))
	    semilogy(xs,s/s[0],'o-',label='singular_value')
	    xlabel('L',size='x-large')
	    ylabel('s_L/s_(L=0)',size='x-large')
	    legend(loc=0)
	    fs=open("singular_value.dat",'w')
	    for i,sv in enumerate(s):
		print>>fs,i,sv
	    fs.close()

	    figure(2)
	    xmesh=linspace(-1,1,U.shape[0])
	    for i in range(Lmax): 
	       plot(xmesh,U[:,i]/U[-1,i],label="L="+str(i))
	    xlabel('x',size='x-large')
	    ylabel('u_L(x)/u_L(1)',size='x-large')
	    legend(loc=0)
	     
	    figure(3)
	    ymesh=linspace(-1,1,V.shape[0])
	    for i in range(Lmax): 
	       plot(ymesh,V[i,:]/V[i,-1],label="L="+str(i))
	    xlabel('y',size='x-large')
	    ylabel('v_L(y)/v_L(1)',size='x-large')
	    legend(loc=0)
	    return (U,s,V)
	        
	    

def safe_exp(x):
	if x<-50.: return 0.0
	return exp(x)
	
def Gl_eva(U,gtau,svd_max):
	Ntau=len(gtau)
	Gl=zeros(svd_max)
	for i in range(svd_max):
	    for j in range(Ntau):
		Gl[i]+=U[j,i]*gtau[j]
	figure(4)
	plot(range(svd_max)[::2],Gl[::2]/Gl[0],'o-',label="even L")
	plot(range(svd_max)[1::2],Gl[1::2]/Gl[0],'s-',label="odd L")
	xlabel('L',size='x-large')
	ylabel('GL',size='x-large')
	legend(loc=0)
	return Gl

def Gtau_fromGl(U,Gl,gtau,tau,Lcut):
	print "svd_lmax=",Lcut,"used to reconstruct G(tau)" 
	Gtau_gl=[]
	Ntau=U.shape[0]
	for i in range(Ntau):
	    sumGtau=0.0
	    for j in range(Lcut):
	        sumGtau+=Gl[j]*U[i,j]
	    Gtau_gl.append(sumGtau)
	Gtau_gl=array(Gtau_gl)
	figure(5)
	semilogy(tau,-gtau,'-',lw=2,label="input")
	semilogy(tau,-Gtau_gl,'--',lw=3,label="IR")
	xlabel(r'$\tau$',size='x-large')
	ylabel(r'$-G(\tau)$',size='x-large')
	legend(loc=0)
	
	
if __name__=="__main__":
	fname="Gtau_qh.dat"
	Gdata=loadtxt(fname).transpose()
	
	Gt=Gtau(Gdata)
	gtau=Gt.gt()
	tau=Gt.tau()
	beta=Gt.beta()
	print "beta=",beta 

	wmax=5;
	lmbda=beta*wmax;ratio=2
	fkgen=fkernel_linearmesh(Gdata)
	kernM=fkgen.fkernel(lmbda,ratio)

	svd_out=svd_de(kernM)
	Lmax=3
	U,s,V=svd_out.svd_eva(Lmax)

	svd_L=20
	Gl=Gl_eva(U,gtau,svd_L)
	
	svd_Lmax=19
	if svd_Lmax>svd_L:
	   print "svd_Lmax has to be smaller than svd_L.."
	else:
	   Gtau_fromGl(U,Gl,gtau,tau,svd_Lmax)
	show()
	
	
	
	
