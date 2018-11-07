#!/usr/bin/env python
from scipy import *
from pylab import *
import glob



if __name__=="__main__":
	dir=glob.glob("mu*")
	
	figure(0)
	for fn in dir:
	    dos=loadtxt(fn+"/dos/dos.out").transpose()
	    plot(dos[0],dos[1],label=fn)
	legend(loc=0)
	xlabel(r'$\omega$',size='x-large')
	ylabel('$DOS$',size='x-large')
	axvline(x=0,c='k',ls='--')
	xlim([-4.5,4.5])

	figure(1)
	for fn in dir:
	    dos=loadtxt(fn+"/dos/dos.out").transpose()
	    plot(dos[0],dos[1],label=fn)
	legend(loc=0)
	xlabel(r'$\omega$',size='x-large')
	ylabel('$DOS$',size='x-large')
	axvline(x=0,c='k',ls='--')
	xlim([-1.5,1.5])

	show()
