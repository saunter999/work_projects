#!/usr/bin/env python
from scipy import *
from pylab import *
from scipy import integrate

def func(x):
	if x<beta/2:
	   return exp(-x)
	else:
	   return -cosh(x)
#	return x 
#	return cosh(x)
#	return exp(x)

def Samplesign(Nitt,warm,measure,tau_init):
	sgn=0.0
	Naver=0.0
	tau=tau_init
	for itt in range(Nitt):
	    taup=rand()*beta
	    P=min([abs(func(taup))/abs(func(tau)),1])
	    if P>rand():
	       tau=taup	

	    if itt>warm and itt%measure==0:
		sgn+=func(tau)/abs(func(tau))
		Naver+=1.0
	return sgn/Naver
		
def Samplenorm(Nitt,warm,measure,Dorder_init):
	N0=0
	N1=0
	Dorder=Dorder_init
	Naver=0.0
	acc1=0.0
	acc2=0.0

	for itt in range(Nitt):
	    if Dorder==0:
		   tau=rand()*beta
		   P=min([beta*abs(func(tau))/Norm,1])
		   if P>rand():#Metropolis
	 	    #accept:
			Dorder=1
			acc1+=1
		    #reject -> do nothing
	    else:
		   P=min([Norm/(beta*abs(func(tau))),1])
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
	Norm=1.0
	print "beta=",beta,"Norm=",Norm
	warm=5000
	Ns=int(1e6)
	Ne=int(8e6)
        Nls=arange(Ns,Ne,Ne/5)	
	res=[]
	measure=100
	t1ls=linspace(0,beta/2,1000)
	t2ls=linspace(beta/2,beta,1000)
	y1ls=array([func(tau) for tau in t1ls]) 
	y2ls=array([func(tau) for tau in t2ls]) 
	exact=integrate.simps(y1ls,t1ls)+integrate.simps(y2ls,t2ls)

	for Nitt in Nls:
	  Dorder_init=0   ##we start at order 0
	  tau_init=rand()*beta
	  N0p,N1p=Samplenorm(Nitt,warm,measure,Dorder_init)
	  sign=Samplesign(Nitt,warm,measure,tau_init) 
	  print "Nitt=",Nitt
	  print "Exact result of the integration is:"
	  print exact
	  print
	  print "Monte carlo sampling is:"
	  print "norm is:",Norm*N1p/N0p,'sign is:',sign
	  print 'result is:',Norm*N1p/N0p*sign
	  res.append(Norm*N1p/N0p*sign)
	  print
	  print "Direct Monte carlo sampling is:"
	  print SimpleMC(Nitt)
	res=array(res)
	plot(Nls,res,'o-')
	axhline(y=exact)
	  
	savefig("General_res.png")
	show()
