#!/usr/bin/env python
# @ Kristjan Haule 2007-2017
from scipy import *
from pylab import *
from scipy import interpolate
import os
import glob
import subprocess

"""
This module runs ctqmc impurity solver for plaquette-DMFT of the Hubbard model.
The executable should exist in directory params['exe']
"""

Uc = 6 # in terms of hopping t
cix_name = 'hubbard_U_'+str(Uc)+'_normal.cix'
print 'cix_name=', cix_name
if not glob.glob(cix_name):
    print """First perform exact diagonalization of the cluster and generate the cix file for this U with the command
     hub -U """+str(Uc)+"""
    """
    sys.exit(0)

# Set mpi_prefix.dat to contain the command for parallel execution, for example "mpirun -n 200"
mpifile = 'mpi_prefix.dat'
MPI = ''
if os.path.isfile(mpifile):
    MPI = open(mpifile, 'r').next().strip()+' '
    
params = {"exe":   ["./ctqmc",          "# Path to executable"],
          "mode":  ["GH",               "# G stands for Green's function sampling, H stands for the Hubbard-I tail"],
          "M" :    [10e6,               "# Number of Monte Carlo steps"],
          "U":     [0,                  "# Should be set to zero because it is already included in cix-file"],
          "mu":    [Uc/2.,              "# Chemical potential"],
          "beta":  [100,                "# Inverse temperature"],
          "nom":   [120,                "# number of Matsubara frequency points to sample"],
          "svd_lmax":[25,                 "# number of SVD functions to project the solution"],
          "tsample":[200,               "# how often to record the measurements" ],
          "aom":   [1,                  "# number of frequency points to determin high frequency tail"],
          "cix":   [cix_name,           "# Input file with atomic state"],
          "Delta": ["Delta.inp",        "# Input bath function hybridization"],
          "GlobalFlip":[1000000,        "# how often to perform global flip"],
          "minF":  [1e-8,               "# anhilation operator is computed only up to this precision."],
          "sderiv":[0.05,               "# the continuity of derivative when merging with the tail"],
	  "PMove": [0.,                 "# we will siwtch off exchange move"],
          }
    
def CreateInputFile(params):
    " Creates input file (PARAMS) for CT-QMC solver"
    f = open('PARAMS', 'w')
    print >> f, '# Input file for continuous time quantum Monte Carlo'
    for p in params:
        print >> f, p, params[p][0], '\t', params[p][1]
    f.close()

    
def AssembleHk(ikx,iky,c,s):
    """This is tight-binding Hamiltonian for 2x2 paquette cellular-DMFT of the Hubbard model
    with nearest neighbor hopping. See PRB 76, 104509 (2007), Eq: 13
    """
    Hk = zeros((4,4), dtype=complex)    
    Hk[0,0] =  -(2.0 + c[ikx] + c[iky])
    Hk[0,1] = s[ikx]*1j
    Hk[0,2] = s[iky]*1j
    Hk[0,3] = 0.0
    Hk[1,1] = c[ikx]-c[iky]
    Hk[1,2] = 0.0
    Hk[1,3] = s[iky]*1j
    Hk[2,2] = -(c[ikx]-c[iky])
    Hk[2,3] = s[ikx]*1j
    Hk[3,3] = 2.0 + c[ikx] + c[iky]
    for i in range(4):
        for j in range(i):
            Hk[i,j] = conjugate(Hk[j,i])
    return Hk

def create_log_mesh(om, nom, ntail_):
    """Creates logarithmic mesh on Matsubara axis
       Takes first istart points from mesh om and the rest of om mesh is replaced by ntail poinst redistribued logarithmically.
       Input:
           om      -- original long mesh
           nom     -- number of points not in the tail
           ntail   -- tail replaced by ntail points only
       Output:
           ind_om  -- index array which conatins index to kept Matsubara points
    """
    istart = min(nom, len(om))
    ntail = min(ntail_, len(om)-istart)
    
    istart = min(nom,len(om))
    ntail = min(ntail, len(om)-istart)

    ind_om=[]
    alpha = log((len(om)-1.)/istart)/(ntail-1.)
    for i in range(istart):
        ind_om.append(i)
    for i in range(ntail):
        t = int(istart*exp(alpha*i)+0.5)
        if (t != ind_om[-1]):
            ind_om.append(t)
    return ind_om

def ComputeGloc(klist,w,Sg,mu,c,s):
    """ Given the set of k-points and self-energies, it computes the local Green's function
    """
    Gloc = zeros((len(w),4,4), dtype=complex)
    Id = identity(4)
    for ikx,iky in klist:
        Hk = AssembleHk(ikx,iky,c,s)
        for iw in range(len(w)):
            G_inv = (w[iw]*1j+mu)*Id - Hk - Sg[iw,:,:] 
            Gloc[iw] += linalg.inv( G_inv )
    Gloc *= 1./len(klist)
    return Gloc

def ComputeDelta(w,Sg,Eimp):
    """ Given the local Green's function, self-energy and impurity levels,
    computes hybridization Delta
    """
    Delta = zeros( (len(w),4,4), dtype=complex)
    Id = identity(4)
    for iw in range(len(w)):
        G_inv = linalg.inv(Gloc[iw,:,:])
        Delta[iw] = w[iw]*1j*Id - Eimp - Sg[iw,:,:] - G_inv
    return Delta

def ReconstructOnEntireMesh(Dx,om,nom,w):
    """ Points in the tail are computed on logarithmic grid. Now we restore
        the tail on all frequency points
    """
    Dspliner = interpolate.UnivariateSpline(w[nom:], Dx[nom:].real, s=0)
    Dsplinei = interpolate.UnivariateSpline(w[nom:], Dx[nom:].imag, s=0)
    Dtail = Dspliner(om[nom:]) + Dsplinei(om[nom:])*1j
    Dt = hstack( (Dx[:nom], Dtail) )
    return Dt

if __name__ == '__main__':
    
    nom =  params['nom'][0]
    mu = params['mu'][0]
    beta = params['beta'][0]
    Nx=24     # mumber of k-points in each direction, in metallic phase you will need to increase this number
    ntail = 30 # number of points in the tail
    EimpCF = [-2,0,0,2]   # crystal fields, given by (\sum_k H_k) of the model
    # Number of DMFT iterations
    Niter = 20

    # Creating parameters file PARAMS for qmc execution
    CreateInputFile(params)

    # impurity levels
    Eimp = zeros((4,4))
    for i in range(4):
        Eimp[i,i] = EimpCF[i]-mu
        
    # create k-mesh
    km = arange(-Nx/2,Nx/2)*2*pi/Nx
    c = cos(km)
    s = sin(km)
    klist = []
    for ikx in range(Nx):
        for iky in range(Nx):
            klist.append( (ikx,iky) )
        
    # current best approximation for the self-energy
    if not glob.glob('Sig.out'): # does not exist
        om = (2*arange(2000)+1.)*pi/beta
        #ss = zeros((len(om),6)) # just set self-energy to zero
        gamma=5e-2 # small broadening
        ss = vstack( (zeros(len(om)), -ones(len(om))*gamma, zeros(len(om)), -ones(len(om))*gamma, zeros(len(om)), -ones(len(om))*gamma) )
        data = vstack( (om, ss) )
    else:
        data = loadtxt('Sig.out').transpose()
    
    for it in range(Niter):
        # self-energy components from file
        om = data[0]
        S00 = data[1]+data[2]*1j
        Sp0 = data[3]+data[4]*1j
        Spp = data[5]+data[6]*1j
        
        # create log-mesh in frequency
        ind_om = create_log_mesh(om,nom,ntail)
        w = zeros(len(ind_om))
        Sg = zeros( (len(ind_om),4,4), dtype=complex)
        for i,iw in enumerate(ind_om):
            w[i] = om[iw]
            # self-energy for plaquette on logarithmic mesh
            Sg[i] = [[S00[iw],0,0,0],[0,Sp0[iw],0,0],[0,0,Sp0[iw],0],[0,0,0,Spp[iw]]]
        
        # Computing the local Green's function
        Gloc = ComputeGloc(klist,w,Sg,mu,c,s)
        # and using the DMFT SCC computing Delta
        Delta = ComputeDelta(w,Sg,Eimp)
        
        # Now interpolate tails on entire mesh
        D00 = ReconstructOnEntireMesh(Delta[:,0,0], om,nom,w)
        Dp0 = ReconstructOnEntireMesh(0.5*(Delta[:,1,1]+Delta[:,2,2]), om,nom,w)
        Dpp = ReconstructOnEntireMesh(Delta[:,3,3], om,nom,w)
        
        datao = vstack( (om,D00.real,D00.imag,Dp0.real,Dp0.imag,Dpp.real,Dpp.imag) )
        savetxt('Delta.inp',datao.transpose())
        
        # Running ctqmc
        print 'Running ---- qmc itt.: ', it, '-----'
        subprocess.call(MPI+params['exe'][0], shell=True,stdout=sys.stdout,stderr=sys.stderr)
    
        # Some copying to store data obtained so far (at each iteration)
        cmd = 'cp Gf.out Gf.out.'+str(it)
        subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr)  # copying Gf
        cmd = 'cp Sig.out Sig.out.'+str(it)
        subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying Sig
        cmd = 'cp ctqmc.log ctqmc.log.'+str(it)
        subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying log file
        cmd = 'cp Delta.inp Delta.inp.'+str(it)
        subprocess.call(cmd, shell=True,stdout=sys.stdout,stderr=sys.stderr) # copying Delta
