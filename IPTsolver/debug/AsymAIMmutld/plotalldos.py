#!/usr/bin/env python
from scipy import *
from pylab import *
import glob
from scipy import integrate



if __name__=="__main__":
	dir=glob.glob("mutld*")
	
	muls=[0.05,-0.05,0.0,-0.1,0.1]
	dosmu=[]
	for i,fn in enumerate(dir):
	    dos=loadtxt(fn+"/dos/dos.out").transpose()
	    plot(dos[0],dos[1],label=fn)
	    print fn
	    dosmu.append(dos)
	legend(loc=0)
	for i in range(len(dir)):
	    dos=dosmu[i]
	    print 'mu=',muls[i]
	    muint=[];dosint=[]
	    for j in range(len(dos[0])):
	        if dos[0][j]<muls[i]:
		   muint.append(dos[0][j])
		   dosint.append(dos[1][j])
	    muint=array(muint)
	    dosint=array(dosint)
	    print integrate.simps(dosint,muint)
	    print
	xlabel(r'$\omega$',size='x-large')
	ylabel('$DOS$',size='x-large')
	xlim([-1.,1.])
	show()
