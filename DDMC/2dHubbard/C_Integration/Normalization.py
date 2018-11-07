#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import integrate

def func(x):
	"""function has to be postive definite in the range [0,beta].
	"""
	return exp(-x)
#	return x 
#	return cosh(x)
#	return exp(x)

def Sample(Nitt,warm,measure,Dorder_init):
	N0=0
	N1=0
	Dorder=Dorder_init
	Naver=0.0
	acc1=0.0
	acc2=0.0

	for itt in range(Nitt):
	    if Dorder==0:
		   tau=rand()*beta
		   P=min([beta*func(tau)/Norm,1])
		   if P>rand():#Metropolis
	 	    #accept:
			Dorder=1
			acc1+=1
		    #reject -> do nothing
	    else:
		   P=min([Norm/(beta*func(tau)),1])
		   if P>rand():#Metropolis
	 	    #accept:
			Dorder=0
			acc2+=1
		    #reject -> do nothing
#	    print acc1,acc2,itt

	    if itt>warm and itt%measure==0:
	        if Dorder==0:
	 	   N0+=1.
	        if Dorder==1:
		   N1+=1.
	        Naver+=1.0
	print "-----------------------------"
	print "number of N0||","number of N1||","sum of these two||"
	print N0,"(",N0/Naver,")",N1,"(",N1/Naver,")",Naver
	print "acceptance rate is:",acc1/Nitt,acc2/Nitt
	return N0/Naver,N1/Naver
		 


def SimpleMC(Nitt):
	sum=0.0
	for itt in range(Nitt):
	     tau=rand()*beta
	     ftau=func(tau)
	     sum+=ftau
	
	return sum/Nitt*beta

if __name__=="__main__":
	"""
	Testing idea of Normalization using the example Int_0^beta d(tau) exp(-tau) 
	"""
	beta=1.5
	Norm=0.1
	print "beta=",beta,"Norm=",Norm
	warm=5000
	Ns=int(1e6)
	Ne=int(8e6)
        Nls=arange(Ns,Ne,Ne/5)	
	res=[]
	measure=100
	tls=linspace(0,beta,1000)
	yls=array([func(tau) for tau in tls]) 
	exact=integrate.simps(yls,tls)

	for Nitt in Nls:
	  Dorder_init=0   ##we start at order 0
	  N0p,N1p=Sample(Nitt,warm,measure,Dorder_init)
	   
	  print "Nitt=",Nitt
	  print "Exact result of the integration is:"
	  print exact
	  print
	  print "Monte carlo sampling is:"
	  print Norm*N1p/N0p
	  res.append(Norm*N1p/N0p)
	  print
	  print "Direct Monte carlo sampling is:"
	  print SimpleMC(Nitt)
	res=array(res)
	plot(Nls,res,'o-')
	axhline(y=exact)
	  
	savefig("res.png")
	show()
