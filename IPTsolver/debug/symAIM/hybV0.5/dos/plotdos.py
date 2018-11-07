#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import integrate



if __name__=="__main__":
	dos=loadtxt("dos.out").transpose()
	plot(dos[0],dos[1])
	xlabel(r'$\omega$',size='x-large')
	ylabel('$DOS$',size='x-large')
	xlim([-3.,3.])
	print "Total spectal weight is:",integrate.simps(dos[1],dos[0])
	show()
