#!/usr/bin/env python
# Author: Kristjan Haule, March 2007-2017
from scipy import * 
from pylab import *
import os,sys,subprocess

##imported class ###
import Giom_to_Gtau
from Giom_to_Gtau import Giom
from Giom_to_Gtau import AsympinverseGiom
from Giom_to_Gtau import Gtau

"""
This module runs ctqmc impurity solver for AIM with flat band dos.
The executable should exist in directory params['exe']
"""

def IPTinput_gen(Nom,beta,ef,mutld):
	"""
	First generate Delta(iom) for the case of flat band dos
	then generate weiss field G_0(iom)
	
	"""
	V,D =0.1,1
	f1 = open("Delta.inp", 'w')
	f2 = open("G0.inp", 'w')
	Delta=[]
	for n in range(Nom):
	    iom = (2.*n+1)*pi/beta
	    Redelta=-V**2/(2.0*D)*log( ((D/2.-mutld)**2+iom**2)/((D/2.+mutld)**2+iom**2))
	    Imdelta=-V**2/D*(arctan((D/2.-mutld)/iom)-arctan((-D/2.-mutld)/iom))
	    print>>f1,iom,Redelta,Imdelta
	    G0= 1./( 1.0j*iom-(ef-mutld)-(Redelta+1.0j*Imdelta) )
	    print>>f2,iom,G0.real,G0.imag
	f1.close()
	f2.close()

if __name__=="__main__":
	beta=50 
	Nom=2000
	ef=0.0;mutld=-0.1
	mu=1.0
	print "beta=",beta
	IPTinput_gen(Nom,beta,ef,mutld)

	print "Reading G0(iom).inp file"
	Gdata=loadtxt("./G0.inp").transpose()
	g0iom=Giom(Gdata)
	G0inv=AsympinverseGiom(Gdata)
	tilte=G0inv.tilte()
	g0tau=Gtau(Gdata)
	np=3.0
	g0tau.Gtau_eva(tilte,np)

	show()

    


        
