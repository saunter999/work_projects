#!/usr/bin/env python
from scipy import *
from pylab import *





if __name__=="__main__":
	"""
	Atomic limit, at half filling we calcuate G_int and G0, and then calcuate self energy
	"""
	Uls=array([0.0,0.5]);beta=50
	muls=Uls/2.
	omls=[];GUom=[];G0om=[];Nom=2000
	for i in range(Nom):
	    om=(2*i+1)*pi/beta
	    omls.append( om )
	    #GUom.append(0.5*(1./(1.0j*om+muls[1])+1./(1.0j*om+muls[1]-Uls[1]/2.)))
	    GUom.append(0.5*(1./(1.0j*om+muls[1])+1./(1.0j*om+muls[1]-Uls[1])))
	    G0om.append(0.5*(1./(1.0j*om+muls[0])+1./(1.0j*om+muls[0]-Uls[0])))
	GUom=array(GUom)
	G0om=array(G0om)
	figure(0)
	plot(omls,GUom.real,label="real"+"U="+str(Uls[1]))
	plot(omls,GUom.imag,label="imag"+"U="+str(Uls[1]) )
	plot(omls,G0om.real,label="real"+"U="+str(Uls[0]))
	plot(omls,G0om.imag,label="imag"+"U="+str(Uls[0]) )
	xlim([0,omls[-1]/20])
	xlabel(r'$i\omega$',size='x-large')
	ylabel(r'$G(i\omega)$',size='x-large')
	legend(loc=0.0)

	figure(1)
	Sigma=[]
	fout=open("Sig.out",'w')
	for i,om in enumerate(omls):
	     sig=1.0/G0om[i]-1.0/GUom[i]
	     print>>fout,om,sig.real,sig.imag 
	     Sigma.append(sig)
	Sigma=array(Sigma)
	fout.close()
	plot(omls,Sigma.real,'o-',lw=2,label="real")
	plot(omls,Sigma.imag,'^-',lw=2,label="imag")
	xlim([0,10])
	xlabel(r'$i\omega$',size='x-large')
	ylabel(r'$\Sigma(i\omega)$',size='x-large')
	legend(loc=0.0)
	
	show()
	
