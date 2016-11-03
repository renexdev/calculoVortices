#!/usr/bin/python
"""
Rene Cejas Bolecek
reneczech@gmail.com / ncejas@cab.cnea.gov.ar
Low Temperatures Laboratory, Centro Atomico Bariloche
Instituto Balseiro, Universidad Nacional de Cuyo
Avenida Bustillo 9500, 8400 Bariloche, Rio Negro, Argentina
MIT License
"""
"Interaction Energy/Force plot script v 5.0"
import numpy as np
import matplotlib.pyplot as plt
from numpy import zeros, sqrt, where, pi, average, arange
import pylab
import matplotlib.ticker
import matplotlib as mpl


################################################################
#CUTOFF energy values greater than threshold_val 
# ***** thdown ---- thup ****
# ----- thup ***** thdown -----
################################################################
def threshold_val(modo_th,val,thr_val_down,thr_val_up):
	if (modo_th == "ext"):
		bools1 = val>thr_val_up
		bools2 = val<thr_val_down
		criterion = bools1+bools2
		
	if (modo_th == "int"):
		bools1 = val<thr_val_up
		bools2 = val>thr_val_down
		criterion = bools1*bools2
	if ((modo_th != "int") and (modo_th != "ext")):
		raise  RuntimeError ("Select correct mode: ext or int")
		
	
	interior_indices, = where(criterion)
	indices_exp = where(criterion)
	num_interior_particles = len(interior_indices)
	if num_interior_particles < 1:
		raise  RuntimeError ("No particles found for which a circle of radius rMax\
		will lie entirely within a square of side length S.  Decrease rMax\
		or increase the size of the square.")
	return (indices_exp)

################################################################
#SETUP SCRIPT

filename = "test_Eint_Fint"
#INPUT FILE
#x_um[i03],y_um[i03],Eint[i03],Eint_E0[i03],fint_x[i03],fint_y[i03],fint_mod.get_length()
# data[:,2] Eint
# data[:,3] Eint/E0
data   = np.loadtxt(filename+'.dat')

var = "E" 
#var: F or E
##############################################################
#mode_trh = "int" or "ext"
mode_trh = "int"
##Lo seleccionado es el rango elegido ( fuera de rango lo pinta de blanco)

thr_in_0=1.18

##5kG  Ethr  Fthr
#@10G:  1.6   4.5/100.0
#@20G:  2.5   9.0/100.0
#@60G:  3.7   44.0/100.0

#thr_in_E =1.6
#thr_in_F =4.5/100.0
#thr_in_E =2.5
#thr_in_F =9.0/100.0
thr_in_E =1.5#3.7
thr_in_F =44.0/100.0

##5kG  Emin  Emax  Fmin  Fmax
#@10G:  0.5   2.3   0   9.0/100.0
#@20G:  1.3   3.5   0   21.0/100.0
#@60G:  2.3   5.1   0   65.0/100.0

################################################################
x = data[:,0]
y = data[:,1]
if(var=="E"):
	z = data[:,3]
	printv='einte0'
	z_index = 3
	thr_in=thr_in_E
	
if(var=="F"):
	z = data[:,6]
	z = z/100
	printv='forcemod'
	z_index = 6
	thr_in=thr_in_F
if(var=="E"):
	#Tick label Controll
	zdelta=0.02#0.001
	
	#zinf=0.5#z.min()
	#zsup=2.3#z.max()
	#zinf=1.3#z.min()
	#zsup=3.5#z.max()
#ESTE ESPECIFICA RANGO DE COLORES
	zinf=0#z.min()
	zsup=2#z.max()
	
	#Color Map Controll
	cinf=zinf#z.min()
	csup= zsup#z.max()
	
if(var=="F"):
	#Tick label Controll
	zdelta=0.02/100#0.001
	zinf=0
		
	#zsup=float(9.0/100.0)
	#zsup=float(21.0/100.0)
	zsup=float(65.0/100.0)
	
	#Color Map Controll
	cinf= zinf
	csup= zsup

print z.min(), z.max()
################################################

indices_exp_out = threshold_val(mode_trh,z,thr_in_0,thr_in)
len_vortices = len(indices_exp_out[0][:])

print "# vortices: ",len(z)
print "# vortices E thr: ",len(indices_exp_out[0][:])

x_um_trh = zeros(len_vortices)
y_um_trh = zeros(len_vortices)
z_thr = zeros(len_vortices)


for i01 in range(len_vortices):
	x_um_trh[i01] = data[indices_exp_out[0][i01],0]
	y_um_trh[i01] = data[indices_exp_out[0][i01],1]
	z_thr[i01] = max(z)
################################################


levels = arange(zinf, zsup, zdelta)
pylab.scatter(x,y,marker='o',c=z,s=35,cmap=pylab.cm.get_cmap('jet', len(levels)-1),norm=mpl.colors.Normalize(vmin=cinf, vmax=csup),linewidths=0) # no 


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

deltap=(zsup-zdelta-(zinf+0.005))/5.0


tick_locs   = [zinf+0.005,zinf+deltap,zinf+2*deltap,zinf+3*deltap,zinf+4*deltap,zsup-zdelta]
cbar = pylab.colorbar( format='%.2f') # draw colorbar
cbar.locator     = matplotlib.ticker.FixedLocator(tick_locs)
cbar.update_ticks()

print x_um_trh
pylab.scatter(x_um_trh,y_um_trh,s=30,marker='o',c='white',linewidths=0,) 

if(var=="E"):
	cbar.set_label(r'$\varepsilon_{int} / \varepsilon_0$', fontsize=30) # _{int} / \varepsilon _{0}
if(var=="F"):
	cbar.set_label(r'$\mid f_i \mid \, (10^{-5}N/m)$', fontsize=30) # _{int} / \varepsilon _{0}

# plot data points.
pylab.xlim(0,max(data[:,0]))
pylab.ylim(0,max(data[:,1]))
pylab.xlabel("$x\, (\mu m)$", fontsize=30)
pylab.ylabel("$y\, (\mu m)$", fontsize=30)
pylab.title('#selected vortices: %d (%d)\n \n'%(len(indices_exp_out[0][:]),len(z)))
pylab.savefig(filename+"_"+printv+'_'+'_thres_'+'%s_'%mode_trh+'%.2f'%thr_in_0+"_"+'%.2f'%thr_in+".png")
pylab.show()


