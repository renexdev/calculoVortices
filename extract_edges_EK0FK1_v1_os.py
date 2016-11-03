#!/usr/bin/python
"""
Rene Cejas Bolecek
reneczech@gmail.com / ncejas@cab.cnea.gov.ar
Low Temperatures Laboratory, Centro Atomico Bariloche
Instituto Balseiro, Universidad Nacional de Cuyo
Avenida Bustillo 9500, 8400 Bariloche, Rio Negro, Argentina
MIT License
"""
"Extract vortices from edge from Energy calculation file"
import pylab
from numpy import zeros, sqrt, where
import numpy as np

def f_a_triang(Bint):	#INPUT: Oe
		y2 = 1.075*sqrt(Phi0/(Bint*1.0))
		return y2
		#OUTPUT: um

##################################################################
Phi0 = 20.7
#[phi0] = Gauss*um^2
#phi0 = 2.07*10^-7Gauss*cm*(10^8um/cm);
Lx_px_in = 540
Ly_px_in = 405

filename = "test_Eint_Fint"

px2um = 0.114791667 #2000 x

Bint_in = 15
##################################################################
dataIn   = np.loadtxt(filename+'.dat')

field_mode = "std" #DEFINITION OF BCS LIMITS

#lattice parameter
atf_px = 11;
afield_px = f_a_triang(Bint_in)/px2um
#how many lattice parameters from edge
param_bcs = 1.01 # dimensionless

a_px = afield_px   # atf_px or afield_px
if(field_mode=="std"):
	
	offset_bcs_px = a_px*param_bcs
	#px units
	Lix_in=(offset_bcs_px)*px2um
	Lfx_in=(Lx_px_in-offset_bcs_px)*px2um
	Liy_in=(offset_bcs_px)*px2um
	Lfy_in=(Ly_px_in-offset_bcs_px)*px2um

###################################################################
def correct_boundaryeffect(x,y,Lix,Lfx,Liy,Lfy):
    bools1 = x>Lix
    bools2 = x<Lfx
    bools3 = y>Liy
    bools4 = y<Lfy
    criterion = bools1*bools2*bools3*bools4
    interior_indices, = where(criterion)
    indices_exp = where(criterion)
    num_interior_particles = len(interior_indices)
    if num_interior_particles < 1:
        raise  RuntimeError ("No vortices found")
	
    return (indices_exp)
###################################################################

x_um = dataIn[:,0]
y_um = dataIn[:,1]


##########################################################################################

indices_exp_out = correct_boundaryeffect(x_um,y_um,Lix_in,Lfx_in,Liy_in,Lfy_in)
len_vortices = len(indices_exp_out[0][:])


x_um_bcs = zeros(len_vortices)
y_um_bcs = zeros(len_vortices)
eint_bcs = zeros(len_vortices)
eine0_bcs= zeros(len_vortices)
fx_bcs= zeros(len_vortices)
fy_bcs= zeros(len_vortices)
fmod_bcs = zeros(len_vortices)
f_write_01 = open(filename+'_BCS_x0_%d'% Lix_in+'_x_%d'% Lfx_in+'_y0_%d'% Liy_in+'_y_%d'% Lfy_in+".dat", "w")

for i04 in range(len_vortices):
	x_um_bcs[i04] = dataIn[indices_exp_out[0][i04],0]
	y_um_bcs[i04] = dataIn[indices_exp_out[0][i04],1]
	eint_bcs[i04] = dataIn[indices_exp_out[0][i04],2]
	eine0_bcs[i04]= dataIn[indices_exp_out[0][i04],3]
	fx_bcs[i04]= dataIn[indices_exp_out[0][i04],4]
	fy_bcs[i04]= dataIn[indices_exp_out[0][i04],5]
	fmod_bcs[i04] = dataIn[indices_exp_out[0][i04],6]
	f_write_01.write("%.2f	%.2f	%.2f	%.2f	%.2f	%.2f\n" % (x_um_bcs[i04],y_um_bcs[i04],eine0_bcs[i04],fx_bcs[i04],fy_bcs[i04],fmod_bcs[i04]))
f_write_01.close()
##########################################################################################

