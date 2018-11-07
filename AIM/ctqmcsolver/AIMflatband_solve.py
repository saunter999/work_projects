#!/usr/bin/env python
# Author: Kristjan Haule, March 2007-2017
from scipy import * 
import os,sys,subprocess
import oneorbcix
from oneorbcix import*

"""
This module runs ctqmc impurity solver for AIM with flat band dos.
The executable should exist in directory params['exe']
"""

Uc,V,D = 3.5,0.2,1
ef,beta=-Uc/2.0,50 
mu=Uc/2.
obcix_gen(ef,Uc)


params = {"exe":   ["./ctqmc",          "# Path to executable"],
          "U":     [0.0,                 "# Coulomb repulsion (F0)"],
          "mu":    [mu,              "# Chemical potential"],
          "beta":  [beta,                "# Inverse temperature"],
          "M" :    [5e6,                "# Number of Monte Carlo steps"],
          "mode":  ["GH",               "# S stands for self-energy sampling, M stands for high frequency moment tail"],
          "cix":   ["one_band.cix",     "# Input file with atomic state"],
          "Delta": ["Delta.inp",        "# Input bath function hybridization"],
          "tsample":[200,               "# how often to record the measurements" ],
          "nom":   [80,                 "# number of Matsubara frequency points to sample"],
        "svd_lmax":[30,                 "# number of SVD functions to project the solution"],
          "aom":   [1,                  "# number of frequency points to determin high frequency tail"],
          "GlobalFlip":[1000000,         "# how often to perform global flip"],
          }



" Creates input file (PARAMS) for CT-QMC solver"
f = open('PARAMS', 'w')
print >> f, '# Input file for continuous time quantum Monte Carlo'
for p in params:
    print >> f, p, params[p][0], '\t', params[p][1]
f.close()

# Preparing input file Delta.inp
f = open(params["Delta"][0], 'w')
for n in range(2000):
    iom = (2.*n+1)*pi/params['beta'][0]
    Redelta=-V**2/(2.0*D)*log( ((D/2.-mu)**2+iom**2)/((D/2.+mu)**2+iom**2))
    Imdelta=-V**2/D*(arctan((D/2.-mu)/iom)-arctan((-D/2.-mu)/iom))
    print>>f,iom,Redelta,Imdelta
f.close()

    




print 'Running ---- qmc ----'
subprocess.call(params['exe'][0], shell=True,stdout=sys.stdout,stderr=sys.stderr)

# Some copying to store data obtained so far (at each iteration)
cmd = 'cp Gf.out Gf.out.'+str(0)
subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying Gf
cmd = 'cp Sig.out Sig.out.'+str(0)
subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying Sig
cmd = 'cp ctqmc.log ctqmc.log.'+str(0)
subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying log file
    
        
