#!/usr/bin/env python
from scipy import *
from pylab import *




if __name__=="__main__":
	GU=loadtxt("/work/qhan/projects/one_site_model/ctqmc_run/mu0.5U/Gf.out").transpose()
	G0=loadtxt("/work/qhan/projects/one_site_model/ctqmc_run/U=0/725694qhan/Gf.out").transpose()
	ref=loadtxt("./ref_selfenergy.dat").transpose()
	omls=GU[0]
	
	Sigma=[]
	for i,om in enumerate(omls):
	     G0om=G0[1]+1.0j*G0[2]
	     GUom=GU[1]+1.0j*GU[2]
	     sig=1.0/G0om[i]-1.0/GUom[i]
	     Sigma.append(sig)
	Sigma=array(Sigma)
	plot(ref[0],ref[1],'o',label="ref")
	plot(omls,Sigma.imag,label="Dyson")
#	xlim([0,20])
	legend(loc=0)
	show()
