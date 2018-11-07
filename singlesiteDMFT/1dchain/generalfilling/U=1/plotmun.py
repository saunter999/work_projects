#!/usr/bin/env python
from scipy import *
from pylab import *



if __name__=="__main__":
	
	data=loadtxt("Mu_n.dat").transpose()
	ref=loadtxt("refU=1_mun.dat").transpose()
	figure(0)
	plot(data[0],data[1],'go-',markersize=9,label="cal")
	plot(ref[0],ref[1],'rs-',markersize=6,label="ref")
	legend(loc=0)

	xlabel(r'$\mu$',size='x-large')
	ylabel('$n$',size='x-large')
#	axvline(x=0,c='k',ls='--')
#	xlim([-4.5,4.5])

	show()
