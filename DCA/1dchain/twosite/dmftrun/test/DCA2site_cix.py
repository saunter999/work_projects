#This file is to generate the cix file for a two impurities Anderson model in which the two impurities don't interact with each directly..
#07/26/2017 Qiang Han
#//////////////////////////////////#
#Labelling for the 16 atomic states of the two impurites. For '_ _|_ _':the first two entries are for down spin and the last two for up spin;among the down spin and up spin config,the first one stores the occupation of the K0 impurity,the second for the Kpi impurity.
#//////////////////////////////////#
def cix_gen(Eimp_K0,Eimp_Kpi,U):
    sts=['0000','1000','0010','0100','0001','1010','0101','1001','0110','1100','0011','1011','1110','0111','1101','1111']
#print len(sts)
    dic={'0000':1,'1000':2,'0010':3,'0100':4,'0001':5,'1010':6,'0101':7,'1001':8,'0110':9,'1100':10,'0011':11,'1011':12,'1110':13,'0111':14,'1101':15,'1111':16}
    occ=[[0.,0.],[1.,0.],[1.,0.],[0.,1.],[0.,1.],[2.,0.],[0.,2.],[1.,1.],[1.,1.],[1.,1.],[1.,1.],[2.,1.],[2.,1.],[1.,2.],[1.,2.],[2.,2.]] ##every entry is [n_K0,n_Kpi].
    bath_num=4

    def num(s):
	"""
	This computes total N.
	"""
        num=0	
        L=list(s)
    	for i in range(len(L)):
	    if L[i]=='1':
	       num=num+1
    	return num

    def spin_z(s):
	"""
	This computes total Sz.
	"""
        spin_z=0
        L=list(s)
        for i in range(len(L)):
	    if i<=1:
	       sign=-1
	    else:
	       sign=1
            spin_z=spin_z+sign*int(L[i])*0.5
        return spin_z

    def index_table(s,bdex):
	"""
	This computes the index table of each impurity state by the impurity creation operator .
	This is used to calculate Z_atom in ctqmc.
	"""
        L=list(s)
        if L[bdex]=='0':
           L[bdex]='1' 
           s=''.join(L)
           index=dic[s]
        else:
           index=0
        return index

    def E_a(s):
	"""
	This computes Eigenenergy of each impurity state.
	"""
    	L=list(s)
    	numK0=0.
   	numKpi=0.
    	dcK0=0.
    	dcKpi=0.
    	E_a=0.
        for i in range(len(L)):
	    if i%2==0 and L[i]=='1':
	       numK0+=1
	    if i%2==1 and L[i]=='1':
	       numKpi+=1
        if numK0==2:
           dcK0=1.
        if numKpi==2:
           dcKpi=1.
        E_a=numK0*Eimp_K0+numKpi*Eimp_Kpi+U*(dcK0+dcKpi)
        return E_a

    def S_t(s):
	"""
	This computes Stotal.
	"""
        L=list(s)
        S_t=0
        S=[0,0]
        for i in range(2):
	    #print i
            if (int(L[i])+int(L[i+2]))==1:
               S[i]=0.5
      	    else:
	       S[i]=0
        for i in range(2):
            S_t=S_t+S[i]
        return S_t

    def dimen(index):
        if index==0:
	   dim=0
        else:
	   dim=1
        return dim

    def mat_elem(s,bdex,index):
        L=list(s)
        sign=0.
        if index==0:
	   mat_elem=''
        else:
           for i in range(bdex):
	       sign=sign+float(L[i]) 
	   mat_elem=(-1.)**sign
        return mat_elem

    f=open("twoimpurty.imp",'w')
    print >>f,"# Cix file for cluster DMFT with CTQMC" #This line is a must-have.
    print >>f,"#",'cluster size,','number of states,','number of baths,','maximum matrix size'
    print >>f,1,16,4,1
    print >>f,"#",'baths,','dimens,','symmetry,','global flip'
    print >>f,0,1,0,0,"# K0 impurity with down spin"  
    print >>f,1,1,1,1,"# Kpi impurity with down spin"
    print >>f,2,1,0,0,"# K0 impurity with up spin"
    print >>f,3,1,1,1,"# Kpi impurity with up spin"
    print >>f,"#",'cluster energies for unique baths,','eps[k]'
    print >>f,Eimp_K0,Eimp_Kpi
    print >>f,'#','\t','N','\t','K','\t','Sz','\t','size','\t','F_K0^{dn,+}','\t','F_Kpi^{dn,+}','\t','F_K0^{up,+}','\t','F_Kpi^{up,+}','\t','Ea','\t','S'
    for i in range(len(sts)):
        print >>f,i+1,'\t',num(sts[i]),'\t',0,'\t',spin_z(sts[i]),'\t',1,'\t',index_table(sts[i],0),'\t\t',index_table(sts[i],1),'\t\t',index_table(sts[i],2),'\t\t',index_table(sts[i],3),'\t\t',E_a(sts[i]),'\t',S_t(sts[i]),'\t',"#",sts[i] #The '\t'or '\t\t' is also a must-have
    print >>f,"#","start-state",'\t',"end-state",'\t',"dim1",'\t',"dim2",'\t',"matrix element"
    for i in range(len(sts)):
        for bdex in range(bath_num):
	    index=index_table(sts[i],bdex)
	    print >>f,i+1,'\t',index,'\t',dimen(i+1),'\t',dimen(index),'\t',mat_elem(sts[i],bdex,index)
    print>>f,"HB2  # Hubbard-I is used to determine high-frequency"
    print>>f,"# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]"
    print>>f,0,0,0,0,0.0
    print>>f,"#","number of operators needed"
    print>>f,0
    f.close()
    f.close()
