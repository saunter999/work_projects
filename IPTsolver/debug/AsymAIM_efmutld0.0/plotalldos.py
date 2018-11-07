#!/usr/bin/env python
from scipy import *
from pylab import *
import glob
from scipy import integrate



if __name__=="__main__":
	dir=glob.glob("ef*")
	
	dosef=[]
	for i,fn in enumerate(dir):
	    dos=loadtxt(fn+"/dos/dos.out").transpose()
	    plot(dos[0],dos[1],label=fn)
	    dosef.append(dos)
	legend(loc=0)

	for i in range(len(dir)):
	    dos=dosef[i]
	    muint=[];dosint=[]
	    for j in range(len(dos[0])):
	        if dos[0][j]<0.0:
		   muint.append(dos[0][j])
		   dosint.append(dos[1][j])
	    muint=array(muint)
	    dosint=array(dosint)
	    print dir[i]
	    print "occ=",integrate.simps(dosint,muint)
	    print
	xlabel(r'$\omega$',size='x-large')
	ylabel('$DOS$',size='x-large')
	xlim([-3.,3.])
	axvline(x=0.1,c='k',ls='--')
	axvline(x=0.2,c='k',ls='--')
	axvline(x=0.5,c='k',ls='--')
	axvline(x=1.0,c='k',ls='--')
	show()
