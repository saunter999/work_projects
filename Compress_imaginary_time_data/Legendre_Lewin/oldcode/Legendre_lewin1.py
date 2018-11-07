#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy.special import eval_legendre as leg
from scipy import integrate

class Gtau:
	"""
	Class Gtau is to  
	"""
	def __init__(self,fname):	
	    self.fn=fname
	    self.data=loadtxt(fname).transpose()
	    print "Input file is:",str(self.fn),"in instantiating Gtau"
	    print

	def beta(self):
	    return self.data[0][-1]

	def tau(self):
	    return self.data[0]

	def gt(self):
	    return self.data[1]

	def Ntau(self):
	    return self.data.shape[1]
class gl:
	def __init__(self,lmax,tau,gtau):
	    self.Lmax=int(lmax)
	    self.tau=tau
	    self.gtau=gtau

	def gl(self): 
	    Gtau=self.gtau
	    tau=self.tau
	    Lmax=self.Lmax
	    Ntau=len(Gtau)
	    beta=tau[-1]
	    Gl=zeros(Lmax)
	    xtau=[]
	    for t in tau:
	        xtau.append(2.0*t/beta-1.)
	    for i in range(Lmax):
	        y=[]
	        for j in range(Ntau):
	            y.append(leg(i,xtau[j])*Gtau[j])
	        y=array(y)
	        Gl[i]=sqrt(2.*i+1)*integrate.simps(y,tau)
	    return Gl

class Gtau_fromGl:
	def __init__(self,Gl,Lcut,tau):
	    self.Gl=Gl
	    self.Lcut=Lcut
	    self.tau=tau

	def Gtau_Gl(self):
	    tau=self.tau
	    Gl=self.Gl
	    Lcut=self.Lcut
	    beta=tau[-1]
	    Ntau=len(tau)
	    Gtau_gl=[]
	    xtau=[]
	    for t in tau:
	        xtau.append(2.0*t/beta-1.)
	    for i in range(Ntau):
	        sumGtau=0.0
		for j in range(Lcut):
		    sumGtau+=sqrt(2.*j+1)/beta*Gl[j]*leg(j,xtau[i])
		Gtau_gl.append(sumGtau)
	    Gtau_gl=array(Gtau_gl)
	    return Gtau_gl
	        

def printclass_attr(classn):
	print "listing only the attributes in the class of"+str(classn)+":"
#	print classn.__dict__
	print dir(classn)
	print

if __name__=="__main__":
	classls=[Gtau,gl,Gtau_fromGl]
	for classn in classls:
	    printclass_attr(classn)
	Lmax=100;Lcut=100
	fname="Gt.dat"

	Gt=Gtau(fname)
	GL=gl(Lmax,Gt.tau(),Gt.gt())
	figure(1)
	plot(range(Lmax),GL.gl())
	Gtau_gl=Gtau_fromGl( GL.gl(), Lcut, Gt.tau() )
	figure(2)
	semilogy(Gt.tau(),-Gt.gt(),'-',lw=2,label="input")
	semilogy(Gt.tau(),-Gtau_gl.Gtau_Gl(),'--',lw=3,label="Legendre")
	
	show()


