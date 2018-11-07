#!/usr/bin/env python
from scipy import *
from pylab import *
from numpy import linalg as LA
import numpy as np


def check_ortho(v1,v2):
	print "orthonormality check"
	print np.dot(v1,v2)


if __name__=="__main__":
	xdim=2;
	ydim=4;
	a = np.random.randn(xdim, ydim)
	print "testing matrix a is:"
	print a
	print
	U,s,V=LA.svd(a)
	print "shpae of U,s,V:"
	print U.shape,s.shape,V.shape
	print

	S=np.zeros((xdim, ydim))
	S[:len(s),:len(s)]=np.diag(s)
	print "singular value matrix is "
	print S
	print
	A_svd=np.dot(U,np.dot(S,V))
	print 'matrix A_svd construsted from svd is'
	print A_svd
	print 
	print "a-A_svd is"
	print a-A_svd
	check_ortho(U[0,:],U[0,:])

