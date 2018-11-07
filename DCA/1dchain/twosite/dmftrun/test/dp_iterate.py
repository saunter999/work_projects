#!/usr/bin/env python
# Author: KH, March 2007
from scipy import * 
import scipy.weave as weave
import io
import os,sys,subprocess
import band as dp_band
#from dp_band import *
import dp         #Importing dp and importing dp.'submodule' are independent
from dp import *
"""
This module runs ctqmc impurity solver for dp model.
The executable shoule exist in directory params['exe']
"""
##################################
#==    dp-Model parameters      ==#
d=2
Udd=3.0
Upp=0.0*Udd
Upd=0.0*Udd
occ=3.
Ep=-2.
Ed=0.
tpp=1.0
tdd=0.2
tpd=1.0
mixr_mu=0.8
Nk=500
Ne=400
dlt=1.0e-2
Uc=0 # Uc is set to be zero for multi-bands
##################################
##################################
#==      MC parameters      ==#
beta=100
N_mats=2000
nom=100
aom=5
M=20e6
Niter=19
#mu=float(loadtxt("mu.inp"))
#Eimp_p=float(loadtxt("Eimp_p.txt")) 
#Eimp_d=float(loadtxt("Eimp_d.txt")) 
##################################

#### Set MPI prefix ##############
fileMPI = 'mpi_prefix.dat'
if (os.path.exists(fileMPI)):
    print "create mpi_prefix."
    mpi_prefix = ' '.join(loadtxt(fileMPI, dtype = str))
else :
    print "mpi_prefix not found"
    mpi_prefix=''
#################################
#########0.Settingband parameters######
dp_band.dp_band.input_params(N_mats,beta,d,Nk,Ne,dlt,Udd,Upp,Upd,occ,Ep,Ed,tpp,tdd,tpd,mixr_mu)

########### 1.kmesh####################    
print "1.calcuating kmesh for irreducible k points"
dp_band.dp_band.klist_gen()
dp_band.dp_band.klist_check()
dp_band.dp_band.kmesh()


###########2.mu_fix####################
print "2.calculating mu for nonint bands"
murng=[-5,5]  #Guess for the right bound for the mu(chemical potential).
mu=dp_band.dp_band.mu_fix(murng)
#f=open("mu.inp",'w')
#print>>f,mu
#f.close()	


########### 3.dos####################    
print "3.calcuating dos for nonint bands"
dp_band.dp_band.dos(mu)

############4.Impurity_level_fix######
print "4.calculating impurity level with double counting correction but without mu"
Eimp=dp_band.dp_band.eimp_fix()


params = {"exe":   [mpi_prefix+" ./ctqmc",       "# Path to executable"],
          "Delta": ["Delta.inp",        "# Input bath function hybridization"],
          "cix":   ["dp.cix",     "# Input file with atomic state"],#cix prepared separtely.
          "U":     [0,                 "# Coulomb repulsion (F0)"],
          "sderiv":[0.05,		"# Maximum derivative mismatch accepted for tail concatenation"],
	  "mu":    [mu,              "# Chemical potential"], #Determined from 1.mu_fix
          "beta":  [beta,                "# Inverse temperature"],
          "M" :    [M,                "# Number of Monte Carlo steps"],
          "nom":   [nom,                 "# number of Matsubara frequency points to sample"],
          "aom":   [aom,                  "# number of frequency points to determin high frequency tail"],
          "tsample":[10,               "# how often to record the measurements" ]}

def CreateInputFile(params):
    " Creates input file (PARAMS) for CT-QMC solver"
    f = open('PARAMS', 'w')
    print >> f, '# Input file for continuous time quantum Monte Carlo'
    for p in params:
        print >> f, p, params[p][0], '\t', params[p][1]
    f.close()
#==============================================================================#
	
def DMFT_SCC(fDelta,mu):  #qh:fDelta is a string specifying the name of file containing hybridization funciton.
    """This subroutine creates Delta.inp from Gf.out for DMFT on 2d square lattice for dp model: Delta_(d/p)=iw-E_imp_(p/d)-Sigma_(d/p)-G_(d/p)^{-1}
    If Gf.out does not exist, it creates Gf.out(Gloc) by doing the k summation of the noninteracting lattice Green's function.
    """
    fileSg = 'Sig.out'
    if (os.path.exists(fileSg)): # If output file exists, start from previous iteration
	print 'DMFT scc looping'
	print "Calculating Glatt_int"
        print "mu",mu,"!!!!!!!!!!!!!!!!!!!!!!!!"
	dp.dpcix_gen(Eimp[0],Eimp[1],Upp,Udd,Upd) 
	dp_band.dp_band.glatt_int(mu)
	print "Calculating Delta_int"
	dp_band.dp_band.delta_int(mu)

    else: # otherwise start from Delta.
	print " Preparing input cix"
	dp.dpcix_gen(Eimp[0],Eimp[1],Upp,Udd,Upd) 
        FileD= 'Delta.inp'
        if(os.path.exists(FileD)):
	  print " Using old Delta"
        else:
	  print 'Starting from non-interacting model'
	  print "Calculating Glatt_nonint"
          print "mu",mu
	  dp_band.dp_band.glatt_nonint(mu)
	  print "Calculating Delta_nonint"
	  dp_band.dp_band.delta_nonint(mu)

def Charge_SCC(it,mu):
    fileGf = 'Gf.out'
    if (os.path.exists(fileGf)):
        murng=[mu-5.0,mu+5.0]
        mu=dp_band.dp_band.charge_scc(it,mu,murng)
    else:
        mu="no need to do Charge_SCC for the initial run" 
    return mu



chem=open("mu.dat","w")	
for it in range(Niter):
    if (it==0):
       mu_mix=mu
       print "Before Charge_SCC"
       print "mu",mu
       params["mu"][0]=mu
       CreateInputFile(params)
       DMFT_SCC(params['Delta'][0],mu)   #qh:hybridization is created and updated when iterated
    else:
       mu_mix=Charge_SCC(it,mu_mix)
       print "After Charge_SCC and mixing"
       print "mu",mu_mix
       params["mu"][0]=mu_mix
       CreateInputFile(params)
       #dp_band.dp_band.dc_correction()
       #Eimp=dp_band.dp_band.eimp_update()
       DMFT_SCC(params['Delta'][0],mu_mix)   #qh:hybridization is created and updated when iterated
    print>>chem,it,mu_mix

    # Running ctqmc
    print 'Running ---- qmc itt.: ', it, '-----'
    #print os.popen(params['exe'][0]).read()
    out, err = subprocess.Popen(params['exe'][0], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    #print err 
    print out

    # Some copying to store data obtained so far (at each iteration)
    cmd = 'cp Gf.out Gf.out.'+str(it) #qh:gf.iteration
    print os.popen(cmd).read() # copying Gf
    cmd = 'cp Sig.out Sig.out.'+str(it)
    print os.popen(cmd).read() # copying Sig
    cmd = 'cp Delta.inp Delta.inp.'+str(it)
    print os.popen(cmd).read() # copying Delta
