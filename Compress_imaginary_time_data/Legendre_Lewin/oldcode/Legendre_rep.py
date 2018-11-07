#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy.special import eval_legendre as leg
from scipy import integrate

def Gl_eva(xtau,tauls,Gtau):
	Lmax=100
#	print len(xtau),len(tauls)
	Ntau=len(Gtau)
	Gl=zeros(Lmax)
	for i in range(Lmax):
	    y=[]
	    for j in range(Ntau):
	        y.append(leg(i,xtau[j])*Gtau[j])
	    y=array(y)
	    Gl[i]=sqrt(2.*i+1)*integrate.simps(y,tauls)
	    #for j in range(Ntau):
	#	Gl[i]+=sqrt(2.*i+1)*Dtau*leg(i,xtau[j])*Gtau[j]
	return Gl

def Gtau_fromGl(xtau,Gl,tauls):
	Lcut=100
	print "LegendrePolynominal_lmax=",Lcut,"used to reconstruct G(tau)" 
	Gtau_gl=[]
	Ntau=len(tauls)
	for i in range(Ntau):
	    sumGtau=0.0
#	    for j in range(0,Lcut,2):
	    for j in range(Lcut):
		sumGtau+=sqrt(2.*j+1)/beta*Gl[j]*leg(j,xtau[i])
#		check+=Gl[j]*sqrt(2.*i+1)/beta
	    Gtau_gl.append(sumGtau)
	Gtau_gl=array(Gtau_gl)
	return Gtau_gl
	    

if __name__=="__main__":
	ftau=loadtxt("Gt.dat").transpose()
#	ftau=loadtxt("Gtau_qh.dat").transpose()
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
	Gl=Gl_eva(xmesh,tauls,Gtau)
	#print "Gl in Legendre basis is "
	#print Gl
	figure(1)
	plot(range(len(Gl)),Gl)
	figure(2)
	Gtau_gl=Gtau_fromGl(xmesh,Gl,tauls)
	semilogy(tauls,-Gtau,'-',lw=2,label="input")
	semilogy(tauls,-Gtau_gl,'--',lw=3,label="Legendre")
	#xlim([32,60])
	#ylim([0.021,0.025])
	legend(loc=0)
	show()
