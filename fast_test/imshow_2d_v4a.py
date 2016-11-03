#!/usr/bin/python
"""
Rene Cejas Bolecek
reneczech@gmail.com / ncejas@cab.cnea.gov.ar
Low Temperatures Laboratory, Centro Atomico Bariloche
Instituto Balseiro, Universidad Nacional de Cuyo
Avenida Bustillo 9500, 8400 Bariloche, Rio Negro, Argentina
MIT License
"""
"Interaction Energy/Force plot script v 4.0"
import numpy as np
import matplotlib.pyplot as plt
from numpy import zeros, sqrt, where, pi, average, arange
import pylab
import matplotlib.ticker

################################################################
#CUTOFF energy values greater than threshold_val 

def threshold_val(val,thr_val):
    bools1 = val>thr_val
    #print bools1, bools2,bools3,bools4
    #input de ejemplo tiene 2276
    criterion = bools1
    interior_indices, = where(criterion)
    indices_exp = where(criterion)
    ##print "indices: ",indices_exp
    num_interior_particles = len(interior_indices)
    #print num_interior_particles, interior_indices[5]
    if num_interior_particles < 1:
        raise  RuntimeError ("No vortices found.")
	
    return (indices_exp)

################################################################
#SETUP SCRIPT
filename = "d05_20G_30um_07_5000_ss44_blur_10_v1_Eint_Fint"
#filename = "test_Eint_Fint"
#INPUT FILE
#x_um[i03],y_um[i03],Eint[i03],Eint_E0[i03],fint_x[i03],fint_y[i03],fint_mod.get_length()
# data[:,2] Eint
# data[:,3] Eint/E0
var = "E" 
#var: F or E

##############################################################
data   = np.loadtxt(filename+'.dat')

x = data[:,0] #um
y = data[:,1] #um
if(var=="E"):
	z = data[:,3]
	printv='einte0'
	z_index = 3
	
if(var=="F"):
	z = data[:,6]
	z = z/100 #proper output units
	printv='forcemod'
	z_index = 6
################################################
#Threshold
ee0 = data[:,3]
f_mod = data[:,6]
thr_in = 2
indices_exp_out = threshold_val(ee0,thr_in)
len_vortices = len(indices_exp_out[0][:])
print len_vortices

x_um_trh = zeros(len_vortices)
y_um_trh = zeros(len_vortices)
z_thr = zeros(len_vortices)


for i01 in range(len_vortices):
	x_um_trh[i01] = data[indices_exp_out[0][i01],0]
	y_um_trh[i01] = data[indices_exp_out[0][i01],1]
	z_thr[i01] = max(z)
################################################

zinf=z.min()
zsup=z.max()
zdelta=0.001
csup= z.max()
cinf= z.min()

levels = arange(zinf, zsup, zdelta)
pylab.scatter(x,y,marker='o',c=z,s=35,cmap=pylab.cm.get_cmap('jet', len(levels)-1),linewidths=0) # no interviene, solo me reescalea el eje c!!
#Option_list
#autumn
#(naranja-amarillo)
#bone
#cool
#copper
#flag
#gray
#hot
#hsv
#jet
#pink
#prism
#spring
#summer
#winter
#spectral

deltap=(csup-zdelta-(cinf+0.005))/5.0

tick_locs   = [cinf+0.005,cinf+deltap,cinf+2*deltap,cinf+3*deltap,cinf+4*deltap,csup-zdelta]
cbar = pylab.colorbar( format='%.2f') # draw colorbar
cbar.locator     = matplotlib.ticker.FixedLocator(tick_locs)
cbar.update_ticks()
pylab.scatter(x_um_trh,y_um_trh,marker='o',c=z_thr,linewidths=0)

if(var=="E"):
	cbar.set_label(r'$\varepsilon_{int} / \varepsilon_0$', fontsize=30) # _{int} / \varepsilon _{0}
if(var=="F"):
	cbar.set_label(r'$\mid f_i \mid \, (10^{-5}N/m)$', fontsize=30) # _{int} / \varepsilon _{0}

# plot data points.
pylab.xlim(0,max(data[:,0]))
pylab.ylim(0,max(data[:,1]))
pylab.xlabel("$x\, (\mu m)$", fontsize=30)
pylab.ylabel("$y\, (\mu m)$", fontsize=30)
pylab.savefig(filename+'_thres_'+printv+".png")
pylab.show()


