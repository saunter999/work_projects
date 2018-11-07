#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import integrate

def realintegrand(tau,om):
	return (1./(tau+2)*exp(1.0j*om*tau)).real
def imagintegrand(tau,om):
	return (1./(tau+2)*exp(1.0j*om*tau)).imag


if __name__=="__main__":
	ns=1000;Nom=20;beta=100
	omls=[]
	Greom=[]
	Gimom=[]
	for i in range(ns,ns+Nom):
	   omls.append((2*i+1)*pi/beta)
#	omls=[0]
	for om in omls:
	    Gre=integrate.quad(realintegrand,0,beta,args=(om,))[0]
	    Gim=integrate.quad(imagintegrand,0,beta,args=(om,))[0]
	    Greom.append(Gre)
	    Gimom.append(Gim)
	plot(omls,Greom,'*-',label='real')
	plot(omls,Gimom,'o-',label='imag')
	legend(loc=0)
	show()
	
