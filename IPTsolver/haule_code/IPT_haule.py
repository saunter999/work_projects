#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import weave
from scipy import integrate
import bisect
import sys, optparse

def Dos(x, t=0.25):
    if (abs(x)<(4*t)): 
        return sqrt((4*t)**2-x**2)/(2*pi*(2*t)**2)
    else:
        return 0.0

def ferm(x):
    if x>100: return 0
    if x<-100: return 1
    return 1./(exp(x)+1.)

def CreateFrequencyIndex(om):
    Nw = len(om)
    iy_m_ix = zeros((len(om),len(om)), dtype=int)
    for ix,x in enumerate(om):
        for iy,y in enumerate(om):
            # We know that mesh is equidistant, i.e.,
            #  x = (ix-Nw/2)*L/(Nw/2)
            # hence z==y-x = ((ix-Nw/2) - (iy-Nw/2))*L/(Nw/2)
            # and   iz = (ix-Nw/2) - (iy-Nw/2) + Nw/2
            iz = (iy-Nw/2) - (ix-Nw/2) + Nw/2
            if iz<0 or iz>=len(om): iz=-1
            iy_m_ix[iy,ix] = iz
    return iy_m_ix

def Hilbert(z, om, DOS):
    D0 = 0.0
    # When imaginary part is small, we treat the integral carefully
    if (abs(z.imag)<0.1):
        D0 = Dos(z.real)
    return integrate.trapz( (DOS-D0)/(z-om) ,om) + D0*(log(z-om[0])-log(z-om[-1]))

def KramarsKronig(om, f):
    dh = float((om[-1]-om[1])/(len(om)-1))  # you need to typecast for weave to know the type
    
    logo = zeros(len(om), dtype=float)
    for i in range(1,len(om)-1):
        logo[i] = log((om[-1]-om[i])/(om[i]-om[0]))
    
    deriv = zeros(len(om), dtype=float)
    deriv[0] = (f[1]-f[0])/dh
    deriv[-1] = (f[-1]-f[-2])/dh
    for i in range(1,len(om)-1):
        deriv[i] = (f[i+1]-f[i-1])/(2*dh)
    
    fr = zeros(len(om), dtype=float)
    ## C++ code
    code="""
       using namespace std;
       for (int i=0; i<om.size(); i++){
          double sum1=0;
          for (int j=0; j<om.size(); j++){
             if (i!=j){
                sum1 += (f(j)-f(i))*dh/(om(j)-om(i));
             }else{
                sum1 += deriv(i)*dh;
             }
          }
          fr(i) = (sum1 + f(i)*logo(i))/pi;
       }
    """
    weave.inline(code, ['om','f','fr','deriv','logo','dh','pi'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return fr

def GetAs(om, G0, frm):
    A0p = zeros(len(om), dtype=float)
    A0m = zeros(len(om), dtype=float)
    ## The C++ code
    code="""
       using namespace std;
       for (int i=0; i<om.size(); i++){
          A0p(i) = -G0(i).imag() * frm(i)/pi;
          A0m(i) = -G0(i).imag() * (1-frm(i))/pi;
       }
    """
    weave.inline(code, ['om','A0m','A0p','G0','frm','pi'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    return (A0p, A0m)

def GetPolarisation(om, A0p, A0m, iy_m_ix):
    dh = float((om[-1]-om[1])/(len(om)-1))  # you need to typecast for weave to know the type
    B1 = zeros(len(om), dtype=float)
    B2 = zeros(len(om), dtype=float)
    ## The C++ code
    code="""
       using namespace std;
       for (int ix=0; ix<om.size(); ix++){
          double sum1=0;
          double sum2=0;
          for (int iy=0; iy<om.size(); iy++){
             int iz = iy_m_ix(iy,ix);
             if (iz>=0){
                sum1 += A0m(iy) * A0p(iz);
                sum2 += A0p(iy) * A0m(iz);
             }
          }
          B1(ix) = sum1*dh;
          B2(ix) = sum2*dh;
       }
    """
    weave.inline(code, ['om','A0p','A0m','iy_m_ix','B1','B2','dh'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    return (B1, B2)

def Sopt(om, A0p, A0m, B1, B2, U, iy_m_ix):
    dh = float((om[-1]-om[1])/(len(om)-1))  # you need to typecast for weave to know the type
    iSg = zeros(len(om), dtype=float)
    ## The C++ code
    code="""
       using namespace std;
       for (int ix=0; ix<om.size(); ix++){
          double sum1=0;
          double sum2=0;
          for (int iy=0; iy<om.size(); iy++){
             int iz = iy_m_ix(iy,ix);
             if (iz>=0){
                sum1 += A0p(iy) * B1(iz);
                sum2 += A0m(iy) * B2(iz);
             }
          }
          iSg(ix) = (sum1+sum2)*dh;
       }
    """
    weave.inline(code, ['om','A0p','A0m','iy_m_ix','B1','B2','dh','iSg'],
                 type_converters=weave.converters.blitz, compiler = 'gcc')
    
    iSg *= -U**2*pi
    return iSg


if __name__ == '__main__':
    """ Computes the solution of the DMFT equations using the
        second order perturbation in U (IPT) for the Hubbard model
        at half-filling.
    """
    usage = """usage: %prog [ options ]"""

    parser = optparse.OptionParser(usage)
    parser.add_option("-L", "--Lw",  type="float", default=4,    help="Frequency cutoff.")
    parser.add_option("-N", "--Nw",  type="int",   default=501,  help="Number of frequency points. Must be odd number")
    parser.add_option("-U", "--U",   type="float", default=4.2,  help="Coulomb repulsion")
    parser.add_option("-T", "--T",   type="float", default=0.01, help="Temperature")
    parser.add_option("-K", "--Nk",  type="int",   default=50,   help="Number of momentum points is Nk x Nk. Should be even number")
    parser.add_option("-g", "--gamma",type="float", default=0.04,  help="Broadening of poles")
    parser.add_option("-i", "--Nitt", type="int",   default=200,  help="Number of self-consistent iterations")
    parser.add_option("-m", "--mix",  type="float", default=0.2,  help="Mixing parameter")
    # Next, parse the arguments
    (options, args) = parser.parse_args()

    mix = options.mix
    
    # creating equidistant frequency mesh
    om = linspace(-options.Lw, options.Lw, options.Nw)
    # Fermi function of the mesh
    frm = array([ferm(x/options.T) for x in om])
    # Index array for the mesh
    iy_m_ix = CreateFrequencyIndex(om)

    # Density of states
    DOS = array([Dos(x) for x in om])

    # Starting self-energy
    Sigma = 0.0
    # Starting G0 corresponding to Sigma=0
    G0 = array([Hilbert(x+0.001j, om, DOS) for x in om])
    old_diff = 1e10
    
    for itt in range(options.Nitt):

        (A0p, A0m) = GetAs(om, G0, frm)
        (B1, B2) = GetPolarisation(om, A0p, A0m, iy_m_ix)
        iSg = Sopt(om, A0p, A0m, B1, B2, options.U, iy_m_ix)
        rSg = KramarsKronig(om,iSg)

        Sigma_new = rSg + iSg*1j

        diff = integrate.trapz(abs(Sigma_new-Sigma),om)
        
        if diff<1e-4: break
        if diff>old_diff:
            mix = mix/1.2
        old_diff = diff
        
        print itt, diff, mix
        
        Sigma = Sigma*(1-mix) + Sigma_new*mix
        
        G = array([Hilbert(om[i]-Sigma[i], om, DOS) for i in range(len(om))])

        G0 = 1/(1/G+Sigma)
        #G0 = 1/(om-(0.5)**2*G)

        for i in range(len(om)):# causality might be violated due to numeric inaccuracy
            if G0[i].imag>0:G0[i] = G0[i].real -1e-10j
                
        if itt% 10 ==0:
            plot(om, -imag(G))

    plot(om, -imag(G), 'k', lw=3)
            
    show()
