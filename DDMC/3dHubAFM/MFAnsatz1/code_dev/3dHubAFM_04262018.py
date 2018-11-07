#!/usr/bin/env python
from scipy import *
from pylab import *

def RBZ_gen(N):
	b1=pi*array([1,1,-1])
	b2=pi*array([1,-1,1])
	b3=pi*array([-1,1,1])
	ratio=linspace(0,1,N,endpoint=False)
	RBZ=[]
	for rx in ratio:
            for ry in ratio:
	        for rz in ratio:
		    k=rx*b1+ry*b2+rz*b3
		    RBZ.append(k)
	return RBZ

def ekdisp(RBZ):
	ek=[]
	for k in RBZ:
	    ek.append( -2.0*t*(cos(k[0])+cos(k[1])+cos(k[2])) )
	return ek
	    
def F_0(h):
	print "Evaluating free energy up at zeroth order for","h=",h
	omega=0.0
	for ek in ekls:
           x=sqrt(ek**2+h**2)-mu
           y=-sqrt(ek**2+h**2)-mu
           if x>0:
              omega+=log(1+exp(-beta*x))
           else :
              omega+=-beta*x+log(1+exp(beta*x))
           if y>0:
              omega+=log(1+exp(-beta*y))
           else :
              omega+=-beta*y+log(1+exp(beta*y))
        omega*=-2.0/beta
        omega*=1./Natom
	return omega

def F_1(h):
	print "Evaluating free energy at first order for",'h=',h
	neff_0up=0.0
	neff_0dn=0.0
	for ek in ekls:
	    neff_0up+=Gckup(ek,h,beta)+GcAup(ek,h,beta)
	    neff_0dn+=Gckup(ek,h,beta)-GcAup(ek,h,beta)
	    
	neff_0up*=1./Ns;neff_0dn*=1./Ns
#	print neff_0up,neff_0dn
	neff_0up=-1.0*neff_0up.real;neff_0dn=-1.0*neff_0dn.real
	print h,neff_0up,neff_0dn,neff_0up+neff_0dn
	omega=U*(neff_0up-0.5+h/U)*(neff_0dn-0.5-h/U)
	omega+=h**2/U
	return omega

def Gqp_ktau(e,tau):
        """
        The Green's function of quasiparticle with energy e in (k,tau) representation
        """
        if e>0 : return -exp(-e*tau)/(1.0+exp(-e*beta))
        else   : return -exp(-e*(tau-beta))/(1.0+exp(e*beta))

def Gckup(ek,h,tau):
	Ek=sqrt(ek**2+h**2)
	return 1./(2.0*Ek*(Ek+ek))*((Ek+ek)**2*Gqp_ktau(Ek-mu,tau)+h**2*Gqp_ktau(-Ek-mu,tau))

def GcAup(ek,h,tau):
	Ek=sqrt(ek**2+h**2)
	return h/(2*Ek)*(Gqp_ktau(Ek-mu,tau)-Gqp_ktau(-Ek-mu,tau))


if __name__=="__main__":
	N=20;Ns=N**3;Natom=Ns*2 ##linear size of supercell comsisting of two atoms
	t,U=1.0,4.0
	beta,mu=100.0,0.0
	hls=linspace(-0.5,0.5,20)
	RBZ=RBZ_gen(N)
	ekls=ekdisp(RBZ)
	betals=linspace(2.0,1.6,1)
	##U=3: beta_c in [2.7,3.3]
	for i,beta in enumerate(betals):
	  figure(i)
	  title('beta='+str(beta)+'U='+str(U))
	  omg0=array([F_0(h) for h in hls])
	  omg1=array([F_1(h) for h in hls])+omg0
#  	  plot(hls,omg0,'*-',label='0th order')
	  plot(hls,omg1,'o-',label='1st order')
#	  ylim([-3,0])
	  legend(loc=0)
	show()
