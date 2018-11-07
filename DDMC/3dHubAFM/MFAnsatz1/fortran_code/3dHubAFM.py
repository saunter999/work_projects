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
        omega=omega*(-2.0/beta)
        omega*=1./Natom
#	omega+=h**2/U
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
	print h,neff_0up-neff_0dn,neff_0up+neff_0dn
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


def Grtau_gen(tauls,h):
	G0up=zeros((Ns,2,Ntau),dtype=complex)
	rls=[];ir=0
	for ix in arange(N):
	    for iy in arange(N):
	      for iz in arange(N):
		rls.append([ix,iy,iz])
		for ik,ek in enumerate(ekls): 
		    G0up[ir,0,:]+=(Gckup(ek,h,tauls)+GcAup(ek,h,tauls))*exp(1.0j*(ix*RBZ[ik][0]+iy*RBZ[ik][1]+iz*RBZ[ik][2]))  ##vector sum in python
		ir+=1
	G0up=G0up/Ns
	return G0up,rls
		    
if __name__=="__main__":
	N=4;Ns=N**3;Natom=Ns*2 ##linear size of supercell comsisting of two atoms
	t,U=1.0,5.0
	beta,mu=1.0,0.0
	Ntau=100
	tauls=linspace(0,beta,Ntau)
	print tauls
        exit()
	hls=linspace(-0.5,0.5,20)
	RBZ=RBZ_gen(N)
	ekls=ekdisp(RBZ)
	print ekls
	
#	for i,beta in enumerate(betals):
#	  figure(i)
	title('beta='+str(beta)+'U='+str(U))
	omg0=[];omg1=[]
	for h in hls:
	    omg0.append(F_0(h))
	    omg1.append(F_1(h))
#	    Grtau,rls=Grtau_gen(tauls,h)
	omg0=array(omg0);omg1=array(omg1)
	print omg0
	omg1+=omg0
	plot(hls,omg0,'*-',label='0th order')
	plot(hls,omg1,'o-',label='1st order')
	legend(loc=0)
	show()
