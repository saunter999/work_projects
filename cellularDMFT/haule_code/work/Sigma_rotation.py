#!/usr/bin/env python
from scipy import *
from numpy import linalg as LA




if __name__=="__main__":
	rot=array([[0.25,0.25,0.25,0.25],[0.25,-0.25,0.25,-0.25],[0.25,0.25,-0.25,-0.25],[0.25,-0.25,-0.25,0.25]])
	print LA.det(rot)
	print LA.inv(rot)
