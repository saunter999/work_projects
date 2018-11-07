#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
from mpl_toolkits.mplot3d import Axes3D


def RBZ_gen(N,d):
	kmesh=linspace(-pi,pi,N,endpoint=False)
	FBZ=[];RBZ=[]
	if d==1:
	   FBZ=kmesh
	   for k in FBZ:
	       if cos(k)>0:
		     RBZ.append(k)
	if d==2:
	  for kx in kmesh:
	      for ky in kmesh:
		      FBZ.append([kx,ky])
	  FBZ=array(FBZ)
	  for k in FBZ:
	       if cos(k[0])+cos([k[1]])>0:
		     RBZ.append(k)
	if d==3:
	  for kx in kmesh:
	      for ky in kmesh:
	          for kz in kmesh:
		      FBZ.append([kx,ky,kz])
	  FBZ=array(FBZ)
	  for k in FBZ:
	       if cos(k[0])+cos([k[1]])+cos(k[2])>0:
		     RBZ.append(k)


	return FBZ,RBZ

if __name__=="__main__":
	N=9;d=3
	FBZ,RBZ=RBZ_gen(N,d)
#	print FBZ,RBZ
	print len(FBZ),len(FBZ)/2,len(RBZ)

	RBZ=array(RBZ).transpose()
	if d==2:	
	  plot(RBZ[0],RBZ[1],'o')
	  xlim([-pi,pi])
	  ylim([-pi,pi])
	if d==3:
	  fig = figure()
	  ax = fig.add_subplot(111, projection='3d')
	  ax.scatter(RBZ[0], RBZ[1],RBZ[2])
	  ax.set_xlim([-pi,pi])
	  ax.set_ylim([-pi,pi])
	  ax.set_zlim([-pi,pi])
	  ax.set_xlabel('X Label')
	  ax.set_ylabel('Y Label')
	  ax.set_zlabel('Z Label')

	show()
