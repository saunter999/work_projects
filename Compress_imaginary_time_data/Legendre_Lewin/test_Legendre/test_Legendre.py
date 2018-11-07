#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy.special import eval_legendre as leg


def check_norm(poly):
	sump=0.0
	N=len(poly)
	for i in range(N):
	    sump+=poly[i]**2
	print sump

if __name__=="__main__":
	NL=6  ## number of Legendre polynominals to be plotted
	Lmin=191;Lmax=196
	Nx=200
	Legp=[]
	x=linspace(-1,1,Nx)
	print "Dx=",x[1]-x[0]
	for n in range(Lmin,Lmax):
	    y=leg(n,x)
	    Legp.append(y)
	    plot(x,y,label='L='+str(n))
	legend(loc=0)
	ylim([-1.2,1.2])
	check_norm(Legp[0])
	show()
