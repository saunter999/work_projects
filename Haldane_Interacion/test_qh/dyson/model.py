import weight
import lattice as lat
from weight import UP,DOWN,IN,OUT
import numpy as np
import math
from logger import *

class BareFactory:
    def __init__(self, map, Lat, Hamiltonian, Anneal):
        self.Lat=Lat
        self.__Map=map
        self.__Model=Hamiltonian["Name"]
        self.__Interaction=np.array(Hamiltonian["Interaction"])
        self.__ExternalField=np.array(Hamiltonian["ExternalField"])
        self.__DeltaField=np.array(Anneal["DeltaField"])
        self.__MaxTauBin=self.__Map.MaxTauBin
        self.__Beta=self.__Map.Beta
        if "Hopping" in Hamiltonian:
            self.__Hopping=np.array(Hamiltonian["Hopping"])
        if "ChemicalPotential" in Hamiltonian:
            self.__Mu=np.array(Hamiltonian["ChemicalPotential"])
        if "Description" in Hamiltonian:
            self.__Description=Hamiltonian["Description"]
        else:
            self.__Description=None

        self.NearestNeighbor=[]
        for i in range(Lat.NSublat):
            self.NearestNeighbor.append([])
            for j in range(Lat.NSublat):
                self.NearestNeighbor[i].append([])
        
        self.NextNearestNeighbor=[]
        for i in range(Lat.NSublat):
            self.NextNearestNeighbor.append([])
            for j in range(Lat.NSublat):
                self.NextNearestNeighbor[i].append([])
        

    def Build(self):
        #self.BareG and self.BareW must be reinitialized at every time Build() is called
        
        self.BareG=weight.Weight("SmoothT", self.__Map, "TwoSpins", "AntiSymmetric", "K", "T")
        self.BareW=weight.Weight("DeltaT", self.__Map, "FourSpins", "Symmetric", "R", "T")
        self.Hoppi=weight.Weight("DeltaT", self.__Map, "TwoSpins", "Symmetric", "R", "T")
       
        LatName=self.Lat.Name
        try:
            getattr(self, self.__Model)(LatName)
        except:
            Assert(False, "Model {0} has not been implemented!".format(self.__Model))
        return (self.BareG,self.BareW,self.Hoppi)

    def DecreaseField(self, Anneal):
        #TODO: what if DeltaField/Interval is not an integer?!
        flag=False
        if abs(self.__DeltaField[0])>1e-5:
            for i in range(len(self.__DeltaField)):
                self.__DeltaField[i] += Anneal["Interval"][i]
                Anneal["DeltaField"][i] += Anneal["Interval"][i] 
            flag=True
        log.info(green("ExternalField decreased to: {0}".format(self.__DeltaField)))
        return flag

    def RevertField(self, Anneal):
        for i in range(len(self.__DeltaField)):
            Anneal["Interval"][i]/=2.0
            self.__DeltaField[i] -= Anneal["Interval"][i]
            Anneal["DeltaField"][i] -= Anneal["Interval"][i]
        log.info(green("ExternalField reverted to: {0}".format(self.__DeltaField)))

    #model defintion
    def DiagCount(self, LatName):
        raise NotImplementedError
    def J1J2(self, LatName):
        raise NotImplementedError
    def Heisenberg(self, LatName):
        Assert(len(self.__Interaction)==1, "Heisenberg model only has one coupling!")
        self.__SpinModel(LatName)
    def Hubbard(self, LatName):
        self.__SpinModel(LatName)
    #calculate bareG and bareW
    def __SpinModel(self, LatName):
        Beta=self.__Map.Beta
        u=self.__Mu
        self.BareG.Data=np.zeros(self.BareG.Shape, dtype=complex)
        self.BareW.Data=np.zeros(self.BareW.Shape, dtype=complex)
        self.Hoppi.Data=np.zeros(self.Hoppi.Shape, dtype=complex)
       
        Interaction=list(self.__Interaction)+[0,0,0,0,0]
        J1,J2,J3=Interaction[0:3]
        Hopping=list(self.__Interaction)+[0,0,0,0,0]
        T1,T2,T3=Interaction[0:3]
    
        J_perturbation=None
        if len(Interaction)>2:
            J_perturbation=Interaction[2:]
        
        #Dimension: 2
        spin=self.__Map.Spin2Index(UP,UP)

        if LatName=="Checkerboard":
        #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            self.NearestNeighbor[A][B]=[(0, 0),(0,Ly-1),(Lx-1,0),(Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0),(0,   1),(1,   0),(   1,   1)]

            self.NextNearestNeighbor[A][A]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]
            self.NextNearestNeighbor[B][B]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0)]

            for i in range(2):
                for j in range(2):
                    #J1 interaction A-->B, B-->A
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J1*SS;
                    #J2 interaction A-->A, B-->B
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[:,i,:,j,self.__Map.CoordiIndex(e)]+= J2*SS;

        
        
        elif LatName=="Haldane":
            #NSublat: 2
            Lx,Ly=self.__Map.L
            A,B=0,1
            Fi = np.pi/2
            self.NearestNeighbor[A][B]=[(0, 0), (Lx-1, 0), (Lx-1,Ly-1)]
            self.NearestNeighbor[B][A]=[(0, 0), (1,0), (1,1)]
            
            self.NextNearestNeighbor[A][A]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0),(1,  1),(Lx-1,  Ly-1)]
            self.NextNearestNeighbor[B][B]=[(0, 1),(1,   0),(0,Ly-1),(Lx-1,   0),(1,  1),(Lx-1,  Ly-1)]
            
            #for su in range(2):
            #  self.BareW.Data[:,su,:,su,0]= J1*SS;
            
            #t1 hopping A-->B, B-->A
            t1=1
            for i in range(2):
                for j in range(2):
                    for e in self.NearestNeighbor[i][j]:
                        self.Hoppi.Data[0,i,0,j,self.__Map.CoordiIndex(e)]+= t1;

            #t2 hopping A-->A, B-->B
            t2=0.0
            for i in range(2):
                for j in range(2):
                    for e in self.NextNearestNeighbor[i][j]:
                      if (i==0):
                        #print e
                        if e==(1,0) or e==(0,1) or e==(Lx-1,Ly-1):
                          self.Hoppi.Data[0,i,0,j,self.__Map.CoordiIndex(e)]+= t2*np.exp(1j*Fi);
                        else:
                          self.Hoppi.Data[0,i,0,j,self.__Map.CoordiIndex(e)]+= t2*np.exp(-1j*Fi);
                      else:
                        if e==(1,0) or e==(0,1) or e==(Lx-1,Ly-1):
                          self.Hoppi.Data[0,i,0,j,self.__Map.CoordiIndex(e)]+= t2*np.exp(-1j*Fi);
                        else:
                          self.Hoppi.Data[0,i,0,j,self.__Map.CoordiIndex(e)]+= t2*np.exp(1j*Fi);

            #J1 interaction A-->B, B-->A  
            J1= 3.0
            for i in range(2):
                for j in range(2):
                    #J1 interaction A-->B, B-->A
                    for e in self.NearestNeighbor[i][j]:
                        self.BareW.Data[0,i,0,j,self.__Map.CoordiIndex(e)]+= J1;
                    #J2 interaction A-->A, B-->B
                    for e in self.NextNearestNeighbor[i][j]:
                        self.BareW.Data[0,i,0,j,self.__Map.CoordiIndex(e)]+= J2;
            

       
        else:
            Assert(False, "Lattice {0} has not been implemented yet!".format(LatName))


        #Bare G
        log.info("set Mu={0}, Hopping={1}, and SmoothT Bare G".format(self.__Mu, self.__Hopping))
        Assert(len(self.__ExternalField)>=self.__Map.NSublat, 
                "expect at least {0} externalfield components!".format(self.__Map.NSublat))
        TauGrid=np.array([self.__Map.IndexToTau(t) for t in range(self.__MaxTauBin)])
        WGrid  =np.array([l for l in range(Lx)])
        
        root2=math.sqrt(2)
        # m is the strength of staggered potential
        m = 0
        mi =np.array([[m,0], 
                      [0,-m]])
        self.Hoppi.FFT("K")        
        for po in range(Lx*Ly):
            self.Hoppi.Data[0,:,0,:,po]+= mi
            Ei,Ui=np.linalg.eig(self.Hoppi.Data[0,:,0,:,po])
            Ui0=Ui.conj()
            Ui1=Ui0.T
            for t in range(self.__MaxTauBin): 
              Gi =np.array([[-np.exp(-(Ei[0]-u)*TauGrid[t])*(1-1/(1+np.exp(Beta*(Ei[0]-u)))),0], 
                            [0,-np.exp(-(Ei[1]-u)*TauGrid[t])*(1-1/(1+np.exp(Beta*(Ei[1]-u))))]])
              self.BareG.Data[0,:,0,:,po,t]=np.dot(np.dot(Ui,Gi),Ui1) 
        #---------------------------------calculate Chern number using hopping matrix-------------- 
        #------G^-1(0,k)=u-e--------  
        shapef=[Lx*Ly,2]
        shapeu=[Lx*Ly]
        fi=np.zeros(shapef, dtype=complex)
        Ux=np.zeros(shapeu, dtype=complex)
        Uy=np.zeros(shapeu, dtype=complex)
        sum0=0
        ui =np.array([[u,0], 
                      [0,u]])
        
        for i in range(Lx*Ly):
          gi = ui-self.Hoppi.Data[0,:,0,:,i]     
          ei,vi=np.linalg.eig(gi)
          #print ei
          #print "-------"
          for m in range(2):
            if (ei[m].real>0):
               #print m
               a=m
          fi[i,:]= vi[:,a]

        for j in range(Ly):
          for i in range(Lx):
            i1 = i+1
            if (i1==Lx):
              i1=0
            j1 = j+1
            if (j1==Ly): 
              j1=0
            for m in range(2):
              Ux[j*Ly+i]=Ux[j*Ly+i]+fi[j*Ly+i,m].conjugate()*fi[j*Ly+i1,m]
            Ux[j*Ly+i]=Ux[j*Ly+i]/abs(Ux[j*Ly+i])
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
          
    def ToDict(self):
        points, lines=self.Lat.GetSitesList(SubLatIn=0)
        interaction=self.__GetBareWList()
        return {"Points":points, "Interaction":interaction, "Lines": lines}

    def Plot(self):
        if self.Lat.Dim==2:
            import matplotlib.pyplot as plt
            Assert(self.Lat.Dim==2, "Plot() only works for two dimensional system for now.")
            color=('r','g','b')
            points, _=self.Lat.GetSitesList()
            for vec, coord, sub in points:
                x,y=vec;
                plt.scatter(x,y,s=100,c=color[sub])
                plt.annotate(
                        str(coord),
                        xy = (x, y), xytext = (15, 10),
                        textcoords = 'offset points', ha = 'right', va = 'bottom',
                        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            for vec, coord, sub in self.__GetBareWList():
                start,end=vec
                x,y=zip(start, end)
                plt.plot(x,y,c=color[sub],lw=2)
            plt.show()

    def __GetBareWList(self):
        Spin=self.__Map.Spin2Index(UP,UP)
        offset=np.array(self.__Map.L)/2-1
        size=self.__Map.Vol*self.__Map.NSublat
        BareWList=[]
        Origin=[0 for e in self.__Map.L]
        for sub in self.__Map.GetAllSublatTuple():
            for coord in self.__Map.GetAllCoordi():
                weight=self.BareW.Data[Spin,sub[IN],Spin,sub[OUT],self.__Map.CoordiIndex(coord)]
                if weight*weight>1.0e-10 and sub[IN]==1:
                    n=self.__Map.LatIndex(coord, sub[OUT])
                    # vec is (real vector of in-site, real vector of out-site) 
                    vec=(self.Lat.GetRealVec(Origin,0, sub[IN], offset), \
                                    self.Lat.GetRealVec(coord,0, sub[OUT], offset))
                    coord=(self.__Map.LatIndex(Origin, sub[IN]), \
                                      self.__Map.LatIndex(coord, sub[OUT]))
                    BareWList.append([vec, coord, sub[IN]])
        return BareWList

