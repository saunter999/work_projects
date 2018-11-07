#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA

def RBZ_gen(N,d):
	kmesh=linspace(-pi,pi,N,endpoint=False)
	FBZ=[];RBZ=[]
	if d==1:
	   FBZ=kmesh
	   for k in FBZ:
	       if k+pi>=-pi and k+pi<pi:
	          if abs(k)>abs(k+pi):
		     RBZ.append(k+pi)
	          else:
		     RBZ.append(k)
	if d==2:
	  for kx in kmesh:
	      for ky in kmesh:
		      FBZ.append([kx,ky])
	  FBZ=array(FBZ)
	  Q1=array([pi,pi])
	  Q2=array([pi,-pi])
	  for k in FBZ:
	      kQ1=k+Q1
	      kQ2=k+Q2
	      if kQ1[0]>=-pi and kQ1[0]<=pi and kQ1[1]>=-pi and kQ1[1]<=pi:
		 if LA.norm(k)>LA.norm(kQ1):
		    RBZ.append(kQ1)
		 else:
		    RBZ.append(k)
	      
	      if kQ2[0]>=-pi and kQ2[0]<=pi and kQ2[1]>=-pi and kQ2[1]<=pi:
		 if LA.norm(k)>LA.norm(kQ2):
		    RBZ.append(kQ2)
		 else:
		    RBZ.append(k)

	return FBZ,RBZ

if __name__=="__main__":
	N=4;d=2
	FBZ,RBZ=RBZ_gen(N,d)
	print FBZ,RBZ
	print len(FBZ),len(RBZ)

	RBZ=array(RBZ).transpose()
	
	plot(RBZ[0],RBZ[1],'o')
	xlim([-pi-1,pi+1])
	ylim([-pi-1,pi+1])
	show()
