#!/usr/bin/python
from vector_class02 import *

import matplotlib.pyplot as plt
import pylab
from numpy import zeros, sqrt, where, pi, average, arange, histogram,log
import numpy as np
from scipy import stats, special


##################################################################
flag_debug = True
flag_print = False
filename = 'bphi100G50um.dat' 
a0 = 16.8287#estimado
param = 0.3 #porcentaje de a0

##################################################################
#data = np.loadtxt(filename+'.dat')
data = np.loadtxt(filename)
f_write_01 = open(filename+"_gt_"+"%.1f"%param+"_a0.dat", "w")

x_px = data[:,0]
y_px = data[:,1]

r_i = Vec2d(0.0,0.0)
r_j = Vec2d(0.0,0.0)
r_ij = Vec2d(0.0,0.0)
cumplen = zeros(len(data))

for i01 in range(len(data)):
	cumplen[i01] = True
if (flag_debug):
	print "len: ",len(data)
for i01 in range(len(data)-1):
	if(flag_print):
		print "i01:",i01
	r_i = Vec2d(x_px[i01],y_px[i01])
	if(cumplen[i01]):
		for i02 in range(i01+1,len(data)-1):
			if(flag_print):
				print "i02:",i02
			r_j = Vec2d(x_px[i02],y_px[i02])
			r_ij = r_i-r_j
			if (r_ij.get_length()<=param*a0):
				cumplen[i02]=False
				if (flag_debug):
					print "no cumple:", r_ij.get_length()," lt ",param*a0,"\n"
			#else:
			#	if (flag_debug):
			#		print "cumple:",r_ij.get_length(),"\n"

for i01 in range(len(data)):
	if(cumplen[i01]):
		f_write_01.write("%d %d\n" % (x_px[i01],y_px[i01]))	
	
f_write_01.close()

