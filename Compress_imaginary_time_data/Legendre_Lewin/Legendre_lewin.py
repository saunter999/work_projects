#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy.special import eval_legendre as leg
from scipy import integrate

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

class gl(Gtau):
	"""
	Class gl is to evaluate the coefficents of G(tau) expanded in terms of Legendre polynominals. 
	Description of its attribute:
	Inheritated attribute from class Gtau
	c.gl_eva(Lmax): Given Lmax--largest order of expansion,it calculate expansion coefficents stored in 1d array of dimension Lmax.
	"""
	def __init__(self,G):
	    Gtau.__init__(self,G)

	def gl_eva(self,Lmax): 
	    Gtau=self.gt()
	    tau=self.tau()
	    Ntau=self.Ntau()
	    beta=self.beta()
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

class Gtau_fromGl(Gtau):
	"""
	Class Gtau_fromGl is to caluate G(tau) with a given set of coefficent of Gl 
	Description of its attribute:
	Inheritated attribute from class Gtau
	Gtau_Gl :Given Lcut--cut off of Lmax and Gl, we calcuate G(tau) 
	"""
	def __init__(self,G):
	    Gtau.__init__(self,G)

	def Gtau_Gl(self,Gl,Lcut):
	    tau=self.tau()
	    beta=self.beta()
	    Ntau=self.Ntau()
	    Gtau_gl=[]
	    xtau=[]
	    Lmax=len(Gl)
	    if Lmax<Lcut:
	       raise Exception('Lcut used to construct Gtau is greater than Lmax!!!')
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
#	"""
#	Note that attributes defined in __init__ part(not in the class of inheritance ) is not counted when list attributes.
#	"""
	print "All attributes(not including inheritead) in the class of"+str(classn)+":"
	print classn.__dict__.keys()
	print "All attributes(including inheritead)  in the class of"+str(classn)+":"
	print dir(classn)
	print

if __name__=="__main__":
	"""
	Reference: PHYSICAL REVIEW B 84, 075145 (2011)
	"""
	classls=[Gtau,gl,Gtau_fromGl]
	for classn in classls:
	    printclass_attr(classn)
	Lmax=50;
	fname="Gt.dat"
	Gdata=loadtxt(fname).transpose()

	Gt=Gtau(Gdata)
	GL=gl(Gdata)
	Glcoef=GL.gl_eva(Lmax)
	figure(1)
	plot(range(Lmax)[::2],Glcoef[::2],'o-',label="even L")
	plot(range(Lmax)[1::2],Glcoef[1::2],'s-',label="odd L")
	xlabel('L',size='x-large')
	ylabel('G_L',size='x-large')
	legend(loc=0)

	figure(2)
	semilogy(Gt.tau(),-Gt.gt(),'-',lw=2,label="input")
#	plot(Gt.tau(),Gt.gt(),'-',lw=2,label="input")
	Gtau_gl=Gtau_fromGl(Gdata)
	Lcutls=[5,10,20]
	for Lcut in Lcutls:
	    Gtau_rec_gl=Gtau_gl.Gtau_Gl(Glcoef,Lcut)
	    semilogy(Gt.tau(),-Gtau_rec_gl,'--',lw=3,label="Lcut="+str(Lcut))
#	    plot(Gt.tau(),Gtau_rec_gl,'--',lw=3,label="Lcut="+str(Lcut))

	xlabel(r'$\tau$',size='x-large')
	ylabel(r'$-G(\tau)$',size='x-large')
	legend(loc=0)
	show()


