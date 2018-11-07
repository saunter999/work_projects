#!/usr/bin/env python
from scipy import *
from scipy import optimize
import scipy.weave as weave
import io
import os,sys,subprocess
import DCA2site_cix
from DCA2site_cix import *

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
	      "cix":   ["twoimpurty.imp",     "# Input file with atomic state"],
	      "U":     [0.0,                 "# Coulomb repulsion (F0)"],
	      "mu":    [mu,              "# Chemical potential"],
	      "beta":  [beta,                "# Inverse temperature"],
	      "M" :    [M,                "# Number of Monte Carlo steps"],
              "mode":  ["GH",               "# S stands for self-energy sampling, M stands for high frequency moment tail"],
              "svd_lmax":[30,                 "# number of SVD functions to project the solution"],
	      "nom":   [2*beta,                 "# number of Matsubara frequency points to sample"],
	      "aom":   [1,                  "# number of frequency points to determin high frequency tail"],
	      "tsample":[300,               "# how often to record the measurements" ],
	      "maxNoise":[1e100,              "# maximum allowed noise is large in simple run"]}
    return params


def CreateInputFile(params):
    " Creates input file (PARAMS) for CT-QMC solver"
    f = open('PARAMS', 'w')
    print >> f, '# Input file for continuous time quantum Monte Carlo'
    for p in params:
        print >> f, p, params[p][0], '\t', params[p][1]
    f.close()

def DMFT_SCC():
    """This subroutine creates Delta.inp for two impurity model.
    If Sig.out exists,Delta=iw+mu-Eimp-Sig-Glatt**(-1)
    If Sig.out does not exist, it creates Delta.inp at the first interation, which corresponds to the non-interacting limit: G_nint=\sum_k 1/(iw+mu-ek) and Delta=iw+mu-Eimp-G_nint**(-1)

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
    gbarK0=zeros(Nmats,dtype=complex)
    gbarKpi=zeros(Nmats,dtype=complex)
    code="""
	# line 75 "2siteDCAhubchain.py"
	using namespace std;
	double efK0=EK0imp;
	double efKpi=EKpiimp;
	for (int i=0;i<masbf.size();i++){
	    for(int j=0;j<Ekls_K0.size();j++){
		gbarK0(i)=gbarK0(i)+1.0/( masbf(i)-Ekls_K0(j) );
		gbarKpi(i)=gbarKpi(i)+1.0/( masbf(i)-Ekls_Kpi(j) );
	    }
	    gbarK0(i)=gbarK0(i)/totk;
	    gbarKpi(i)=gbarKpi(i)/totk;
	    deltaK0(i)=masbf(i)-efK0-1.0/gbarK0(i);
	    deltaKpi(i)=masbf(i)-efKpi-1.0/gbarKpi(i);
	}
	"""
    weave.inline(code,['EK0imp','EKpiimp','masbf','Ekls_K0','Ekls_Kpi','gbarK0','gbarKpi','totk','deltaK0','deltaKpi'],type_converters=weave.converters.blitz,force=1,compiler='gcc')
    for ind,iom in enumerate(masbf):
         print>>f,iom.imag,deltaK0[ind].real,deltaK0[ind].imag,deltaKpi[ind].real,deltaKpi[ind].imag
    f.close()

def Delta_update():
    ##writing new Delta into Delta.inp file
    f=open("Delta.inp",'w')
    gbarK0=zeros(Nmats,dtype=complex)
    gbarKpi=zeros(Nmats,dtype=complex)
    code="""
	# line 101 "2siteDCAhubchain.py"
	using namespace std;
	double efK0=EK0imp;
	double efKpi=EKpiimp;
	for (int i=0;i<masbf.size();i++){
	    for(int j=0;j<Ekls_K0.size();j++){
		gbarK0(i)=gbarK0(i)+1.0/( masbf(i)+mu-Ekls_K0(j)-SigitK0(i) );
		gbarKpi(i)=gbarKpi(i)+1.0/( masbf(i)+mu-Ekls_Kpi(j)-SigitKpi(i));
	    }
	    gbarK0(i)=gbarK0(i)/totk;
	    gbarKpi(i)=gbarKpi(i)/totk;
	    deltaK0(i)=masbf(i)+mu-efK0-1.0/gbarK0(i)-SigitK0(i);
	    deltaKpi(i)=masbf(i)+mu-efKpi-1.0/gbarKpi(i)-SigitKpi(i);
	}
	"""
    weave.inline(code,['EK0imp','EKpiimp','masbf','Ekls_K0','Ekls_Kpi','gbarK0','gbarKpi','totk','deltaK0','deltaKpi','SigitK0','SigitKpi','mu'],type_converters=weave.converters.blitz,force=1,compiler='gcc')
    for ind,iom in enumerate(masbf):
         print>>f,iom.imag,deltaK0[ind].real,deltaK0[ind].imag,deltaKpi[ind].real,deltaKpi[ind].imag
    f.close()

def Gloc_lattice():
    f=open("Glatt_loc.out",'w')
    glatt=zeros(Nmats,dtype=complex)
    code="""
	# line 125 "2siteDCAhubchain.py"
	using namespace std;
	for (int i=0;i<masbf.size();i++){
	    for(int j=0;j<Ekls_K0.size();j++){
		glatt(i)=glatt(i)+1.0/( masbf(i)+mu-Ekls_K0(j)-SigitK0(i) );
		glatt(i)=glatt(i)+1.0/( masbf(i)+mu-Ekls_Kpi(j)-SigitKpi(i));
	    }
	    glatt(i)=glatt(i)/(2.0*totk);
	}
	"""
    weave.inline(code,['masbf','Ekls_K0','Ekls_Kpi','glatt','totk','SigitK0','SigitKpi','mu'],type_converters=weave.converters.blitz,force=1,compiler='gcc')
    for ind,iom in enumerate(masbf):
         print>>f,iom.imag,glatt[ind].real,glatt[ind].imag
    f.close()

def Read_sigout():
    fileSg='Sig.out'
    data=loadtxt(fileSg)
    sig1=zeros(Nmats,dtype=complex)
    sig2=zeros(Nmats,dtype=complex)
    for i in range(Nmats):
        sig1[i]=data[i][1]+1.0j*data[i][2]
        sig2[i]=data[i][3]+1.0j*data[i][4]
    return (sig1,sig2)

#def mu_update(mu_t):
#	Ne=Nmu_compute(mu_t)
#	return Ne-Nocc

#def Nmu_compute(mu_t):
#	Ne=zeros(1)
#	gnontk=0.0j
#        gintk=0.0j
#	code="""
#	#line 157 "1dchainHubiterate.py"
#	using namespace std;
#	#include <math.h> 
#	for (int i=0;i<Ekls.size();i++){
#	      Ne(0)=Ne(0)+1./( exp(beta*(Ekls(i)-mu_t))+1.0);
#	      for(int j=0;j<masbf.size();j++){
#		 gnontk=1.0/( masbf(j)+mu_t-Ekls(i) );
# 	         gintk=1.0/( masbf(j)+mu_t-Ekls(i)-Sigit(j) );
#		 Ne(0)=Ne(0)+2.0/beta* real(gintk-gnontk); 
#	      }
#        }
#	Ne(0)=Ne(0)*2.0/totk;  // 2.0 is due to spin degeneracy.
#	"""
#        weave.inline(code,['Ne','gnontk','gintk','masbf','Ekls','Sigit','mu_t','totk','beta'],type_converters=weave.converters.blitz,compiler='gcc')
#	return Ne[0]


if __name__=="__main__":
    #-------model parameter-----#
    Nocc=1.0;
    t,Uc=1.0,3.0
    Nk=2000 ## 2*Nk is the number of k points in each momentum patch 
    totk=2.0*Nk
    mu=Uc/2.
    mixr=0.9

    #-------MC parameters------#
    beta=100.0;Nmats=2000;
    M=9e6;Niter=19

    mpi_prefix=setMPIprefix()
    params=setparams(mpi_prefix)
    CreateInputFile(params)

    masbf=zeros(Nmats,dtype=complex)
    for i in range(Nmats):
	masbf[i]=1.0j*(2.0*i+1.)*pi/beta
    deltaK0=zeros(Nmats,dtype=complex)
    deltaKpi=zeros(Nmats,dtype=complex)
    SigitK0=zeros(Nmats,dtype=complex)
    SigitKpi=zeros(Nmats,dtype=complex)

    ##1.generate k mesh in each momentum patch
    K0mesh=array(list(linspace(0,0.25,Nk))+list(linspace(0.75,1.0,Nk)))
    Kpimesh=linspace(0.25,0.75,2*Nk)

    ##2.Mapping onto two-impurity AIM model##
    EK0imp=0.0
    EKpiimp=0.0
    Ekls_K0=[]
    Ekls_Kpi=[]
    for k in K0mesh:
	EK0imp+=-2.0*t*cos(k*2.0*pi)
	Ekls_K0.append( -2.0*t*cos(k*2.0*pi) )
    for k in Kpimesh:
	EKpiimp+=-2.0*t*cos(k*2.0*pi)
	Ekls_Kpi.append( -2.0*t*cos(k*2.0*pi) )

    Ekls_K0=array(Ekls_K0);Ekls_Kpi=array(Ekls_Kpi)
    EK0imp=EK0imp/totk
    EKpiimp=EKpiimp/totk
    cix_gen(EK0imp,EKpiimp,Uc)


    fmu=open("mu.dat",'w')
    print>>fmu,"#iteration,mu"
    for it in range(Niter):
	if it==0: print "0th iteration"
	if(it>0):
#	    print it,"th iteration"
            (SigitK0,SigitKpi)=Read_sigout()
	    print it,'Update mu using Glatt(iom)'
	    mu_old=mu
	  #  Ne=Nmu_compute(mu_old)
	   # diff=abs(Ne-Nocc) 
	   # if diff>0.01:
	    #   mu=optimize.brentq(mu_update,mu_old-Uc*diff,mu_old+Uc*diff)
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
   ##calculting Gloc of the lattice using the converged self energy in different patch
    (SigitK0,SigitKpi)=Read_sigout() 
    Gloc_lattice()	
