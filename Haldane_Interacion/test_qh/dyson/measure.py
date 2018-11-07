import numpy as np
import calculator as calc
import lattice as lat
from weight import UP,DOWN,IN,OUT
from logger import *
import os, sys, weight
PI=np.pi

class Observable:
    def __init__(self, Map, lat):
        self.__History={}
        self.__Lat=lat
        self.__Map=Map

    def Append(self, Name, Value):
        if self.__History.has_key(Name) is False:
            self.__History[Name]=[]
        self.__History[Name].append(Value)

    def Measure(self, Chi, Determinate, G, NN):
        self.Append("1-JP", Determinate.min())
        Factor=self.__Map.Beta/self.__Map.MaxTauBin
        if self.__Lat.Name in ["Square", "Cubic", "Checkerboard", "3DCheckerboard", "ValenceBond"]:
            Chi.FFT("R", "W")
            UnifK=(0.0,)*self.__Map.Dim
            #print UnifK
            StagK=(PI,)*self.__Map.Dim
            #print StagK
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, [UnifK,StagK],"Real")
            self.Append("UnifChi", ChiK[0])
            self.Append("StagChi", ChiK[1])
            self.Append("StagMag^2 density", ChiK[1]/self.__Map.Vol/self.__Map.NSublat/self.__Map.Beta)

        elif self.__Lat.Name in ["Pyrochlore"]:
            Chi.FFT("R", "W")
            KList=[(0.0,0.0,0.0), (4*PI, 2*PI ,0)] #High symmetry point with strongest order
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, KList,"Real")
            self.Append("UnifChi", ChiK[0])
            self.Append("Chi_X(4Pi,2Pi,0)", ChiK[1])

        elif self.__Lat.Name in ["Kagome","Haldane"]:
            G.FFT("K","T")
            Lx,Ly=self.__Map.L
            MaxTauBin=self.__Map.MaxTauBin
            TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(MaxTauBin)])
            fillt0=0
            fillt1=0
            for po in range(Lx*Ly):
              fillt0+=G.Data[0,:,0,:,po,0].trace()
              fillt1+=G.Data[0,:,0,:,po,1].trace()
            print "the filling factor is:------------->"  
            print (fillt0-(fillt1-fillt0)/(TauGrid[1]-TauGrid[0])*TauGrid[0])/(2*Lx*Ly)+1   
            dt=TauGrid[1]-TauGrid[0]    
            Beta=self.__Map.Beta
            wn=np.pi/Beta
            shapeg=[2,2]
            shapef=[Lx*Ly,2]
            shapeu=[Lx*Ly]
            fi=np.zeros(shapef, dtype=complex)
            Ux=np.zeros(shapeu, dtype=complex)
            Uy=np.zeros(shapeu, dtype=complex)
            sum0=0
          #---------------------------------calculate Chern number using green function-------------- 
            for i in range(Lx*Ly):
              Gi=np.zeros(shapeg, dtype=complex)
              for t in range(MaxTauBin):
                Gi=Gi+G.Data[0,:,0,:,i,t]*np.exp(1j*TauGrid[t]*wn)*dt
              gi= np.linalg.inv(Gi)
              ei,vi=np.linalg.eig(gi)        
              a=-1
              for m in range(2):
                if (ei[m].real>0):
                  a=m
              if (a==-1):
                 fi[i,:]=np.array([0,0])
              else:
                 fi[i,:]= vi[:,a]
              
            for j in range(Ly):
                for i in range(Lx):
                    i1 = i+1
                    if (i1==Lx):
                       i1=0
                    j1 = j+1
                    if (j1==Lx): 
                       j1=0
                    for m in range(2):
                      Ux[j*Lx+i]=Ux[j*Lx+i]+fi[j*Lx+i,m].conjugate()*fi[j*Lx+i1,m]
                    Ux[j*Lx+i]=Ux[j*Lx+i]/abs(Ux[j*Lx+i])
      
                    for m in range(2):
                      Uy[j*Lx+i]=Uy[j*Lx+i]+fi[j*Lx+i,m].conjugate()*fi[j1*Lx+i,m]
                    Uy[j*Lx+i]=Uy[j*Lx+i]/abs(Uy[j*Lx+i])
            
        
            for j in range(Ly):
                for i in range(Lx):
                    i1 = i+1
                    if (i1==Lx):
                        i1=0
                    j1 = j+1
                    if (j1==Lx): 
                        j1=0
                    mi=1j*(np.log(Ux[j*Lx+i])+np.log(Uy[j*Lx+i1])-np.log(Ux[j1*Lx+i])-np.log(Uy[j*Lx+i]))
                    if mi.real>3.14159 or mi.real<-3.14159:
                         print i+j*Ly, mi.real
                    sum0 = sum0 + 1j*np.log(Ux[j*Lx+i]*Uy[j*Lx+i1]/(Ux[j1*Lx+i]*Uy[j*Lx+i]))
            print "--------"
            print sum0
            Chi.FFT("R", "W")
            KList=[(0.0,0.0),] #High symmetry point with strongest order
            _, ChiK=self.__Lat.FourierTransformation(Chi.Data[0,:,0,:,:,0]*Factor, KList,"Real")
            self.Append("UnifChi", ChiK[0])
        else:
            Assert(False, "model not implemented!")

        if self.__Lat.Name in ["ValenceBond"]:
            Chi.FFT("R","W")
            Neighbor1=self.__Map.CoordiIndex((0,0))
            Neighbor2=self.__Map.CoordiIndex((0,self.__Lat.L[1]-1))
            self.Append("VBS", (Chi.Data[0,0,0,1,Neighbor1,0]-Chi.Data[0,0,0,1,Neighbor2,0])*Factor)

        Chi.FFT("R","W")
        energy=0j
        for i in range(self.__Lat.NSublat):
            for j in range(self.__Lat.NSublat):
                for l in NN[i][j]:
                    energy+=Chi.Data[0,i,0,j,self.__Map.CoordiIndex(l),0]/self.__Lat.NSublat
                    #print self.__Map.CoordiIndex(l)
        self.Append("Energy", energy/self.__Map.MaxTauBin)

        G.FFT("R","T")
        for i in range(self.__Map.NSublat):
            self.Append("<Sz_{0}>".format(i), 0.5*(G.Data[UP,i,UP,i,0,-1]-G.Data[DOWN,i,DOWN,i,0,-1]))
            #print G.Data[UP,i,UP,i,0,-1]
        infostr="Latest measurement:\n"
        for key in sorted(self.__History.keys()):
            infostr+="{0}={1}\n".format(key, self.__History[key][-1])
        log.info(infostr)

    def Load(self, FileName):
        try:
            Dict=IO.LoadDict(FileName)
        except:
            Dict=self.__History  #use old History if we fail to load
        finally:
            self.__History=Dict

    def Save(self, FileName):
        #TODO: add error bar estimation
        IO.SaveDict(FileName, "w", self.__History)
