def obcix_gen(Eimp,U):
    # This file is using binary representation for the atomic states to write a cix file for the oneband Hubbard model
    #Four atomic states are representated by "_ _" where the first entry store the occupation of the down spin:0 means no;1 means one down spin;for the up spin it is similar.
    #===========================================#
    #/////////////////////////////////////////////////////////////////#
    #"sts" stores the foure atomic states:'00'-the vacuum;'10'-one down spin;
    #'01'-one up spin;'11'two opposite spins.                     
    #"dic" stores the index of the foure states.                   
    #"num_bath":number of baths.
    #////////////////////////////////////////////////////////////////#
    sts=['00','10','01','11']
    dic={'00':1,'10':2,'01':3,'11':4}
    num_bath=2
    #===========================================#
    def num(s):
	"""
	This function counts the total number of electrons of the atomic state:
	For '1' in the states in "sts",the num is added by 1.
	"""
	num=0
	f=list(s)
	for i in range(len(f)):
	    if f[i]=='1':
	       num=num+1
	return num
    def spin_z(s):
	"""
	This function computes the total spin in z direction for the atomic state
	"""
	sz=0
	num=0
	f=list(s)
	for i in range(len(f)):
	    if f[i]=='1':
	       sz=sz+(-1)**(i-1)*0.5  #for the '1' in the first entry,the spin in z direction is added by -0.5;for the '1' in the second entry,the spin in z direction is added by 0.5.
	return sz
    def bath(s,dic,bdex):
	"""
	This function computes the result of two baths after actiong on the four atomic states.
	bdex is the index for the two baths:0:bath of down spin;1:bath of up spin.
	"""
	index=0
	f=list(s)
	if f[bdex]=='0':
	    f[bdex]='1'    #after acted upon by the creation operator of the bath of bdex,the number is changed into 1.
	    s=''.join(f)
	    index=dic[s]
	else:
	    index=0
	return index
    def E_a(s,ind):
	E_atom=0.
	L=list(s)
	for i in range(len(L)):
	    E_atom+=Eimp*float(L[i])  
	if ind==3:
	    E_atom+=U          #Coulomb U is included here for one-band model.
	return E_atom
    def S_t(s):
	"""
	This function computes the total spin of a given atomic state
	"""
	st=0
	if s=='00'or s=='11':
	    st=0
	else:
	    st=0.5
	return st
    def dim(index):
	"""
	This function computes the dimension of the start state and end state
	"""
	dimen=0
	if index==0:
	    dimen=0
	else:
	    dimen=1
	return dimen
    def mat_elem(s,index,bdex):
	"""
	This function computes the matrix element of the <end-state|\Psi^{bath,+}|start-state>
	"""
	if index==0:
	   matrix=''    #If the end state is 0,no matrix element.
	else:
	   f=list(s)
	   if f[0]=='1'and bdex==1: #If the down spin is occupied and a creation operator of up spin is acting upon it,we have a minus sign due to fermion statistics.
	      matrix=-1.
	   else:
	      matrix=1.
	return matrix
	    
    f=open("one_band.cix",'w')
    print >>f,"#","Cix file for cluster DMFT with CTQMC"
    print >>f,'#','cluster size,','number of state,','number of baths,','maximum matrix size'
    print >>f,1,4,2,1
    print >>f,'#','baths,','dimens,','symmetry,','global flip'
    print >>f,0,1,0,0
    print >>f,1,1,0,0
    print >>f,'#','cluster energies for unique baths,','eps[k]'
    print >>f,Eimp,0
    print >>f,'#','\t','N','\t','K','\t','Sz','\t','size','\t','F^{dn,+}','\t','F^{up,+}','\t','Ea','\t','S'
    for i in range(len(sts)):
    #    print >>f,i+1,'',num(sts[i]),'',0,'',spin_z(sts[i]),'',1,'',bath(sts[i],dic,0),'',bath(sts[i],dic,1),'',0,'',S_t(sts[i])
	print >>f,i+1,'\t',num(sts[i]),'\t',0,'\t',spin_z(sts[i]),'\t',1,'\t',bath(sts[i],dic,0),'\t\t',bath(sts[i],dic,1),'\t\t',E_a(sts[i],i),'\t',S_t(sts[i])
    print >>f,'#','start-state','\t','end-state','\t','dim1','\t','dim2','\t','matrix elements'
    for i in range(len(sts)):
	for bdex in range(num_bath):
	    index=bath(sts[i],dic,bdex)
	    #print >>f,dic[sts[i]],'',index,'',dim(index),'',dim(index),'',mat_elem(sts[i],index,bdex)
	    print >>f,dic[sts[i]],'\t',index,'\t',dim(index),'\t',dim(index),'\t',mat_elem(sts[i],index,bdex)
    print>>f,"HB2  # Hubbard-I is used to determine high-frequency"
    print>>f,"# UCoulomb : (m1,s1) (m2,s2) (m3,s2) (m4,s1)  Uc[m1,m2,m3,m4]"
    print>>f,0,0,0,0,0.0
    print>>f,"#","number of operators needed"
    print>>f,0
    f.close()
