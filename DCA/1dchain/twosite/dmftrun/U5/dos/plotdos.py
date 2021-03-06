#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import integrate
import glob



if __name__=="__main__":
	"""
	chemical potential has been shifted to 0 in "dos.out" generated by maxent
	"""
	dir=glob.glob("dos.out*")
	for idx,fn in enumerate(dir):
	    figure(idx)
	    dos=loadtxt(fn).transpose()
	    plot(dos[0],dos[1],label="orb"+str(idx))
	    xlabel(r'$\omega$',size='x-large')
	    ylabel('$DOS$',size='x-large')
	    axvline(x=0,c='k',ls='--')
	    xlim([-9.,9.])
	    legend(loc=0)
	    print "Total spectal weight(without spin) is:",integrate.simps(dos[1],dos[0])

	    muint=[];dosint=[]
	    for i in range(len(dos[0])):
		if dos[0][i]<0.0:
		  muint.append(dos[0][i])
		  dosint.append(dos[1][i])
	    muint=array(muint)
	    dosint=array(dosint)
	    print "(n_up+n_dn) is", 2*integrate.simps(dosint,muint) , "in orb "+str(idx) 
	    print
	show()
