#!/usr/bin/env python
from scipy import *
from pylab import *


if __name__=="__main__":
	data=loadtxt("TN_phasediagram.txt").transpose()
	plot(data[0],data[1],'o-')
	xlabel("U")
	ylabel("T_N")
	show()
