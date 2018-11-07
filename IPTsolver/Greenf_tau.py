#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import integrate

class Giom:
	"""
	Class Giom is encapsulate G(iom) in a single class object
	Description of its attribute: 
	c.iom() : 1d array of imaginary frequency data iom
	c.Gre() : 1d array of real part G(iom)
	c.Gim() : 1d array of imag part G(iom)
	c.beta(): inverse temperature beta
	"""
	def __init__(self,Gdata):
	    self.data=Gdata
	def iom(self):
	    return self.data[0]
	def Gre(self):
	    return self.data[1]
	def Gim(self):
	    return self.data[2]
	def beta(self):
	    return pi/self.data[0][0]

class AsympinverseGiom(Giom): 
	"""
	Class AsympinverseGiom is to display the behavior of "G^{-1}(iom)-iom" and extract "tilte" from the tail of "G^{-1}(iom)-iom" 
	Description of its attribute: 
	Inheritated attributes from Giom
	c.tilte() : Return vaule--tilte, which is set as the real part of the tail of "G^{-1}(iom)-iom".
	            Not as a return value--plotting "G^{-1}(iom)-iom" vs "omega". 
	"""
	def __init__(self,G):
	    Giom.__init__(self,G)

	def tilte(self):
	    iom=self.iom()
	    gre=self.Gre()
	    gim=self.Gim()
	    Giom=gre+1.0j*gim
	    ginv=[]
	    for idx,g in enumerate(Giom):
		ginv.append(1./g-1.0j*iom[idx])
	    ginv=array(ginv)
	    figure(1)
	    semilogx(iom,ginv.real,"^-",label="Real")
	    semilogx(iom,ginv.imag,"o-",label="Imag")
	    axhline(y=0,c='k',ls='--')
	    xlabel(r'$i\omega$',size='x-large')
	    ylabel(r'$1/G(i\omega)-i\omega$',size='x-large')
	    legend(loc=0)
	    return ginv[-1].real

class Gtau(Giom):
	"""
	Class Gtau is to calculate G(tau) from G(iom) using modified summation formula.
	Description of the attributes:
	Ihheritated attributes from Giom
	c.Gtau_eva() : Not as a return value--Evaluate G(tau) for a mesh of tau
	"""
	def __init__(self,G):
	    Giom.__init__(self,G)

	def Gtau_eva(self,tilte,np):
	    iom=self.iom()
	    gre=self.Gre()
	    gim=self.Gim()
	    beta=self.beta()
	    Giom=gre+1.0j*gim
	    Niom=len(iom)
	    print "beta=",beta
	    Ntau=int(np*beta)
	    tau=[]
	    lmbda=1.1
	    for i in reversed(range(Ntau)):
	        tau.append(lmbda**(-i)*beta)
	        tau.append(-lmbda**(-i)*beta)
	    tau.sort()
	    print "number of tau points in [-beta,beta]=",len(tau)
	    gtau=[]
	    for t in tau:
		sumG=0.0
		for i in range(Niom): 
#	            sumG+=gre[i]*cos(iom[i]*t)+gim[i]*sin(iom[i]*t)    ##formula without subtraction--this leads to discontinuity at tau=0 and tau=beta
		    tiltG=Giom[i]-1./(1.0j*iom[i]-tilte)
		    sumG+=tiltG.real*cos(iom[i]*t)+tiltG.imag*sin(iom[i]*t)
		sumG=sumG*2./beta
		if t>0:
		   gtau.append( sumG-exp(-t*tilte)/(exp(-beta*tilte)+1.) )
	        else:
		   gtau.append( sumG+exp(-t*tilte)/(exp(beta*tilte)+1.) )
	    gtau=array(gtau)
	    print "occ=",gtau[int(len(tau)/2)-1]
	    figure(2)
	    plot(tau,gtau,'o')
	    #xlim([-beta/2.0,beta/2.0])
	    xlabel(r'$\tau$',size='x-large')
	    ylabel(r'$G(\tau)$',size='x-large')
	    axhline(y=0,c='k',ls='--')
	    axvline(x=0,c='k',ls='--')

	    fout=open("Gtau_qh.dat",'w')
	    print "Generating file 'Gtau_qh.dat'..."
	    print>>fout,"#tau,G(tau)"
	    for i in range(len(tau)):
	          print>>fout,tau[i],gtau[i]
	    fout.close()

def printclass_attr(classn):
	print "All attributes(not including inheritead) in the class of"+str(classn)+":"
	print classn.__dict__.keys()
	print "All attributes(including inheritead)  in the class of"+str(classn)+":"
	print dir(classn)
	print

def Giom_fromGtau(gt):
	ntau=gt.shape[1]
	beta=gt[0][-1]
	Nom=2000
	omls=[]
	Gom=[]
	for i in range(Nom):
	    omls.append( (2*i+1)*pi/beta)

	for om in omls:
	  y=[]
	  for i in range(ntau):
	      y.append( gt[1][i]*exp(1.0j*om*gt[0][i]) )
	  Gom.append(integrate.simps(y,gt[0]))
	Gom=array(Gom)
	figure(3)
	plot(omls,Gom.imag,label='imag')
	plot(omls,Gom.real,label='real')
	xlabel(r'$i\omega$',size='x-large')
	ylabel(r'$G(i\omega)$',size='x-large')
	axhline(y=0,c='k',ls='--')
	axvline(x=1./gt[0][1],c='g',ls='--')
	legend(loc=0)
	fout=open("Gw_qh.dat",'w')
	print "Generating file 'Gw_qh.dat'..."
	print>>fout,"#omega,ReG(i\omega),ImG(i\omega)"
	for i in range(Nom):
	    print>>fout,omls[i],Gom[i].real,Gom[i].imag
	fout.close()
	 

if __name__=="__main__":

	classls=[Giom,AsympinverseGiom,Gtau]
	for classn in classls:
	    printclass_attr(classn)

	print "Reading Gf.out file"
	Gdata=loadtxt("./G0.inp").transpose()
	giom=Giom(Gdata)
	Ginv=AsympinverseGiom(Gdata)
	tilte=Ginv.tilte()

	gtau=Gtau(Gdata)
	np=2.0
	gtau.Gtau_eva(tilte,np) ##np is nunmber of points of tau in unit of beta, i.e # of points is np*beta

##	Gtau_qh=loadtxt("Gtau_qh.dat").transpose()		
##	Giom_fromGtau(Gtau_qh)
	show()
