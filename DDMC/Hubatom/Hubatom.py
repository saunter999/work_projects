#!/usr/bin/env python
from scipy import *
from pylab import *
def F_exact(beta):
	return -U/2.-1/beta*(log(2.)+log(1+exp(-beta*U/2.)))


if __name__=="__main__":
	U=10;betals=linspace(1,100,30)
	F=array([F_exact(beta) for beta in betals])
	plot(betals,F)
	axhline(y=-U/2,c='r')
	ylim([-U/2-1,-U/2+1])
	xlabel('beta')
	ylabel('F')
	show()
