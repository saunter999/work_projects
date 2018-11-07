#!/usr/bin/env python
from scipy import *
from scipy import optimize
import scipy.weave as weave
import io
import os,sys,subprocess

def setMPIprefix():
    fileMPI = 'mpi_prefix.dat'
    if (os.path.exists(fileMPI)):
	print "create mpi_prefix."
	mpi_prefix = ' '.join(loadtxt(fileMPI, dtype = str))
    else :
	print "mpi_prefix not found"
	mpi_prefix=''
    return mpi_prefix

def setparams(mpi_prefix):
    """
	mu-dependent
    """
    params = {"exe":   [mpi_prefix+" ./ctqmc",       "# Path to executable"],
	      "Delta": ["Delta.inp",        "# Input bath function hybridization"],
	      "cix":   ["singleband.imp",     "# Input file with atomic state"],
	      "U":     [Uc,                 "# Coulomb repulsion (F0)"],
	      "mu":    [mu,              "# Chemical potential"],
	      "beta":  [beta,                "# Inverse temperature"],
	      "M" :    [M,                "# Number of Monte Carlo steps"],
              "mode":  ["SM",               "# S stands for self-energy sampling, M stands for high frequency moment tail"],
              "svd_lmax":[30,                 "# number of SVD functions to project the solution"],
	      "nom":   [2*beta,                 "# number of Matsubara frequency points to sample"],
	      "aom":   [1,                  "# number of frequency points to determin high frequency tail"],
	      "tsample":[300,               "# how often to record the measurements" ],
	      "maxNoise":[1e100,              "# maximum allowed noise is large in simple run"]}
    return params

def cixgen():
  icix="""# Cix file for cluster DMFT with CTQMC
  # cluster_size, number of states, number of baths, maximum matrix size
  1 4 2 1
  # baths, dimension, symmetry, global flip
  0       1 0 0
  1       1 0 0
  # cluster energies for unique baths, eps[k]
  0 0
  #   N   K   Sz size F^{+,dn}, F^{+,up}, Ea  S
  1   0   0    0   1   2         3        0   0
  2   1   0 -0.5   1   0         4        0   0.5
  3   1   0  0.5   1   4         0        0   0.5
  4   2   0    0   1   0         0        0   0
  # matrix elements
  1  2  1  1    1    # start-state,end-state, dim1, dim2, <2|F^{+,dn}|1>
  1  3  1  1    1    # start-state,end-state, dim1, dim2, <3|F^{+,up}|1>
  2  0  0  0
  2  4  1  1   -1    # start-state,end-state, dim1, dim2, <4|F^{+,up}|2>
  3  4  1  1    1
  3  0  0  0
  4  0  0  0
  4  0  0  0
  HB2                # Hubbard-I is used to determine high-frequency
  # UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]
  0 0 0 0 0.0
  # number of operators needed
  0
  """
  return icix

def CreateInputFile(params):
    " Creates input file (PARAMS) for CT-QMC solver"
    f = open('PARAMS', 'w')
    print >> f, '# Input file for continuous time quantum Monte Carlo'
    for p in params:
        print >> f, p, params[p][0], '\t', params[p][1]
    f.close()

def DMFT_SCC():
    """This subroutine creates Delta.inp for ctqmc on 1d chain.
    If Sig.out exists,Delta=iw+mu-Eimp-Sig-Glatt**(-1)
    If Sig.out does not exist, it creates Delta.inp at the first interation, which corresponds to the non-interacting limit: G_nint=\sum_k 1/(iw+mu-ek) and Delta=iw+mu-Eimp-G_nint**(-1)
    In the latter case it also creates the impurity cix file, which contains information about the atomic states.

    mu--dependent
    """
    fileSg = 'Sig.out'
    if (os.path.exists(fileSg)): # If Sig.out file exists, start from previous iteration
        print "Updated mu is=",mu
        print "DMFT scc loop~"
	print "Calculating new Delta.inp"
	Delta_update()
    else: # otherwise start from Delta.
        FileD= 'Delta.inp'
        if(os.path.exists(FileD)):
	  print " Using old Delta"
        else:
	  print 'Starting from non-interacting limit'
	  print "Calculating Delta_nonint"
	  Delta_start()

def Delta_start():
    ##writing starting Delta into Delta.inp file(we take mu=0 and sigma=0,we dont need to use mu for a given n..)
    f=open("Delta.inp",'w')
    glatt=zeros(Nmats,dtype=complex)
    code="""
	# line 103 "1dchainHubiterate.py"
	using namespace std;
	double ef=Eimp;
	for (int i=0;i<masbf.size();i++){
	    for(int j=0;j<Ekls.size();j++){
		glatt(i)=glatt(i)+1.0/( masbf(i)-Ekls(j) );
	    }
	    glatt(i)=glatt(i)/totk;
	    delta(i)=masbf(i)-ef-1.0/glatt(i);
	}
	"""
    weave.inline(code,['masbf','Ekls','glatt','delta','Eimp','totk'],type_converters=weave.converters.blitz,force=1,compiler='gcc')
    for ind,iom in enumerate(masbf):
         print>>f,iom.imag,delta[ind].real,delta[ind].imag
    f.close()

def Delta_update():
    ##writing new Delta into Delta.inp file
    f=open("Delta.inp",'w')
    glatt=zeros(Nmats,dtype=complex)
    code="""
	# line 124 "1dchainHubiterate.py"
	using namespace std;
	double ef=Eimp;
	for (int i=0;i<masbf.size();i++){
	    for(int j=0;j<Ekls.size();j++){
		glatt(i)=glatt(i)+1.0/( masbf(i)+mu-Ekls(j)-Sigit(i) );
	    }
	    glatt(i)=glatt(i)/totk;
	    delta(i)=masbf(i)+mu-ef-1.0/glatt(i)-Sigit(i);
	}
	"""
    weave.inline(code,['masbf','Ekls','glatt','delta','Sigit','mu','Eimp','totk'],type_converters=weave.converters.blitz,force=1,compiler='gcc')
    for ind,iom in enumerate(masbf):
         print>>f,iom.imag,delta[ind].real,delta[ind].imag
    f.close()

def Read_sigout():
    fileSg='Sig.out'
    data=loadtxt(fileSg)
    sigma=zeros(Nmats,dtype=complex)
    for i in range(Nmats):
        sigma[i]=data[i][1]+1.0j*data[i][2]
    return sigma

def mu_update(mu_t):
	Ne=Nmu_compute(mu_t)
	return Ne-Nocc

def Nmu_compute(mu_t):
	Ne=zeros(1)
	gnontk=0.0j
        gintk=0.0j
	code="""
	#line 157 "1dchainHubiterate.py"
	using namespace std;
	#include <math.h> 
	for (int i=0;i<Ekls.size();i++){
	      Ne(0)=Ne(0)+1./( exp(beta*(Ekls(i)-mu_t))+1.0);
	      for(int j=0;j<masbf.size();j++){
		 gnontk=1.0/( masbf(j)+mu_t-Ekls(i) );
 	         gintk=1.0/( masbf(j)+mu_t-Ekls(i)-Sigit(j) );
		 Ne(0)=Ne(0)+2.0/beta* real(gintk-gnontk); 
	      }
        }
	Ne(0)=Ne(0)*2.0/totk;  // 2.0 is due to spin degeneracy.
	"""
        weave.inline(code,['Ne','gnontk','gintk','masbf','Ekls','Sigit','mu_t','totk','beta'],type_converters=weave.converters.blitz,compiler='gcc')
	return Ne[0]


if __name__=="__main__":
    #-------model parameter-----#
    Nocc=0.9;
    t,Uc=1.0,7.0
    Nk=2000
    mixr=0.9 ##mu mixing ratio
    totk=(Nk-0.0)**1
    mu=Uc/2. ##initial guess for mu,better choice???

    #-------MC parameters------#
    beta=100.0;Nmats=2000;
    M=5e6;Niter=19

    mpi_prefix=setMPIprefix()
    params=setparams(mpi_prefix)
    CreateInputFile(params)
    icix=cixgen()
    f = open(params['cix'][0], 'w')
    print >> f, icix
    f.close()

    masbf=zeros(Nmats,dtype=complex)
    for i in range(Nmats):
	masbf[i]=1.0j*(2.0*i+1.)*pi/beta
    delta=zeros(Nmats,dtype=complex)
    Sigit=zeros(Nmats,dtype=complex)
    
    ## 1.setting up for noninteracting band dispersion##
    Ekls=[]
    kmesh=linspace(-pi,pi,Nk,endpoint=False) ##important to excluding the endpoint
    for kx in kmesh:
	Ekls.append(-2.0*t*cos(kx))
    Ekls=array(Ekls) 

    
    ##2.Mapping onto AIM model##
    Eimp=0.0
    for ek in Ekls:
	Eimp+=ek
    Eimp=Eimp/totk


    fmu=open("mu.dat",'w')
    print>>fmu,"#iteration,mu"
    for it in range(Niter):
	if it==0: print "0th iteration"
	if(it>0):
	    print it,"th iteration"
            Sigit=Read_sigout()
	    print it,'Update mu using Glatt(iom)'
	    mu_old=mu
	    Ne=Nmu_compute(mu_old)
	    diff=abs(Ne-Nocc) 
	    if diff>0.01:
	       mu=optimize.brentq(mu_update,mu_old-Uc*diff,mu_old+Uc*diff)
	    mu=mu*mixr+mu_old*(1.0-mixr)
	    print>>fmu,it,mu
  	    params=setparams(mpi_prefix)
  	    CreateInputFile(params)
	DMFT_SCC()

	# Running ctqmc
	print 'Running ---- qmc itt.: ', it, '-----'
	#print os.popen(params['exe'][0]).read()

	out, err = subprocess.Popen(params['exe'][0], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
	#subprocess.call(params['exe'][0],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE
	print err 
	print out

	# Some copying to store data obtained so far (at each iteration)
	cmd = 'cp Gf.out Gf.out.'+str(it)
	print os.popen(cmd).read() # copying Gf
	cmd = 'cp Sig.out Sig.out.'+str(it)
	print os.popen(cmd).read() # copying Sig
	cmd = 'cp Delta.inp Delta.inp.'+str(it)
	print os.popen(cmd).read() # copying Delta
    fmu.close()
