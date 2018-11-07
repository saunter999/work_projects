#!/usr/bin/env python
from scipy import *
from pylab import *




if __name__=="__main__":
	ref=loadtxt("./ref_selfenergy.dat").transpose()
	data=loadtxt("/work/qhan/projects/one_site_model/ctqmc_run/mu0.5U/Sig.out").transpose()
##	data=loadtxt("/work/qhan/projects/one_site_model/Sig.out").transpose()
	plot(ref[0],ref[1],'o',label="ref")
	plot(data[0],data[2],label="data")
	xlim([0,20])
	legend(loc=0)
	show()
