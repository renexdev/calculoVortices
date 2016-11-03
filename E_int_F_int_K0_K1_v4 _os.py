#!/usr/bin/python
"""
Rene Cejas Bolecek
reneczech@gmail.com / ncejas@cab.cnea.gov.ar
Low Temperatures Laboratory, Centro Atomico Bariloche
Instituto Balseiro, Universidad Nacional de Cuyo
Avenida Bustillo 9500, 8400 Bariloche, Rio Negro, Argentina
MIT License
"""
"Interaction Energy calculation v 4.0"
from vector_class02_os import *
import matplotlib.pyplot as plt
import pylab
from numpy import zeros, sqrt, where, pi, average, arange, histogram,log
import numpy as np
from scipy import special 
##################################################################
Phi0 = 20.7
#[phi0] = Gauss*um^2
#phi0 = 2.07*10^-7Gauss*cm*(10^8um/cm);
##################################################################

def f_lambda(t,lambda0): #INPUT: dimensionless, um
		y1 = lambda0/sqrt(1-(t*t*t*t))
		return y1
		#OUTPUT: um
def f_a_triang(Bint):	#INPUT: G
		y2 = 1.075*sqrt(Phi0/(Bint*1.0))
		return y2
		#OUTPUT: um
def f_Eint_0(lambdaT):#INPUT: um
		y3 = Phi0*Phi0/(16*pi*pi*lambdaT*lambdaT)
		#Gcm^2=dynacm^2/esu=Oecm^2 
		#G^2cm^4= OeGcm^2= erg / cm^3 . cm^2	= erg/cm
		#OUTPUT(SI): 
		#OUTPUT(cgs): (1e-8)erg/cm=(1e-12)erg/um
		return y3
def f_Eint_log(r,lambdaT):	#INPUT: um, um 
		y4 =Phi0*Phi0/(8*pi*pi*lambdaT*lambdaT)*(log(lambdaT/(r*1.0)) + 0.12)
		return y4
		#Gcm^2=dynacm^2/esu=Oecm^2 
		#G^2cm^4= OeGcm^2= erg / cm^3 . cm^2	= erg/cm
		#OUTPUT(SI): 
		#OUTPUT(cgs): (1e-8)erg/cm=(1e-12)erg/um
def f_Eint_K(r,lambdaT):	#INPUT: um, um 
		y5=Phi0*Phi0/(8*pi*pi*lambdaT*lambdaT)*(special.k0((r*1.0)/lambdaT))
		return y5
		#Gcm^2=dynacm^2/esu=Oecm^2 
		#G^2cm^4= OeGcm^2= erg / cm^3 . cm^2	= erg/cm
		#OUTPUT(SI): 
		#OUTPUT(cgs): (1e-8)erg/cm=(1e-12)erg/um
def f_fint_K1(r_i,r_j,lambdaT,E0):	
		rij_2d = Vec2d(0,0)
		rij_n_2d = Vec2d(0,0)
		y6_2d = Vec2d(0,0)
		rij_2d = r_i-r_j
		rij_n_2d =rij_2d.normalized()
		y6_2d=2*E0/lambdaT*(special.k1(rij_2d.get_length()/lambdaT))*rij_n_2d
		#print y6_2d
		#(10-8 erg / cm ) /um = (10-8 dyne) / (10-4cm) = dyne/cm * 10^-4
		# dyne/cm * 10^-4
		# dyne/cm * 10^-4 = 10^-4 dyne/cm . (10^-5 N/dyne) . (cm/10^-2m) = 10 ^-7 N/m
		#OUTPUT (SI):(1e-7)N/m
		#OUTPUT (cgs): (1e-4)dyn/cm = (1e-8)dyn/um
		return y6_2d	

###################################################################
#Inputs
MAT = "Bi2212"
filename="test"
px2um = 0.114791667 #2000 x
pic_mag = 2000
#FACTOR FOR THE RESIZED IMAGE: YOU MUST CONSERVE ASPECT RATIO	
#=645*(82,65/645)*1500/MAG/540 calibrated for 1500 
#MAG= 4000 -> 30.99/540=0.05738
#229.583333333/MAG
#229.583333333/value
lambda0 = 0.18#(def:0.18)\bibitem{keessuper33200} C.J. van der Beek and M. Konczykowski, R.J. Drost and P.H. Kes, N.Chikumoto, S. Bouffard. Physica C Superconductivity {\b332}, 178 (2000).
atf_px = 10	#PX 10.3
#extrapolated from irreversibility line
Tc = 88.7	#[Tc]=K		
#Hall meas, from fitted curve @9.6G 
T_irr=86.1	#[T_irr]=K
T_melt=T_irr	#[T_melt]=K
H_appl = 15	#[H_appl]=Oe
B_int = 12.3 #from SEM img	#[B_int]=G

###################################################################
data   = np.loadtxt(filename+'.dat')
#[data] = pixel Units

#print px2um,px2um_math 
x_um = data[:,0]*px2um
y_um = data[:,1]*px2um
atf_um = atf_px*px2um # Average lattice parameter calculated from T.F
		
#\lambda_{ab}(0)_um + Tc(affected by irrad)_K + Tim_K + par_n (cut-off radius for calculation)

lambda0_in = lambda0
#[lambda0]=um
par_n = 10 #(default 10)
#cut-off radius for calculation: interaction energy between vortices is negligible when distance between vortices are 
#more than 10 lattice parameter lenght 
rMax_in= par_n*atf_um  #cut-off radius for AV. LATTICE PARAMETER
rMax_trian_in=f_a_triang(B_int)*par_n #cut-off radius for IDEAL LATTICE PARAMETER
print "rMax(um) =",rMax_in, "rMax_triang(um)  =",rMax_trian_in
#outer diameter 

Tc_in = Tc
#[Tc]=K
T_in = T_irr #
#[T_freez]=K
##################################################################

t_in = T_in/(Tc_in*1.0)

#lambdaT_in=f_lambda(t_in,lambda0_in) #2fluid model
#lambdaT_in=1.2 #um #ad hoc value
lambdaT_in=lambda0 #um #experimental lambda
print "Using lambda: ", lambdaT_in, "um"


E00= f_Eint_0(lambdaT_in)
Eint_redttiang =f_Eint_K(rMax_trian_in/par_n,lambdaT_in)
print "E_int_red_perfecta",Eint_redttiang,"  E_int_red_perfecta/e0",Eint_redttiang/E00
#TEST MATH VERSION
#lambdaT_in=lambda0_in


print "$\lambda_0$=",lambda0_in,"[um]\n","$T_c$ =", Tc_in,"[K]\n","$T_{frez}$ =",T_in,"[K]\n","$t_{f}$ =",t_in,"\n","$\lambda(T)$ =",lambdaT_in,"[um]\n","$E_0$ =",f_Eint_0(lambdaT_in),"[(1e-8)erg/cm]"
f_write_01 = open(filename+"_Param_Eint"+".dat", "w")
f_write_01.write("MAT=%s\n#Totalvotices=%d\nlambda0[um]=%.2f\nTc[K]=%.2f\nT[K]=%.2f\nt=%.5f\nlambda(T,Tc)[um]=%.5f\nEnergy_vortex_line(lambda)[(1e-8)erg/cm]=%.5f\nEint_trianglatt(B,lambda)[(1e-8)erg/cm]=%.5f\nEint_trianglatt/E0=%.5f\nHappl[G]=%.2f\nParameter=%i\nrMAX(TF)[um]=%.2f\nrMAX(TRIANG_LAT)[um]=%.2f" % (MAT,len(data),lambda0_in,Tc_in,T_in,t_in,lambdaT_in,f_Eint_0(lambdaT_in),Eint_redttiang,Eint_redttiang/E00,H_appl,par_n,rMax_in,rMax_trian_in))
f_write_01.close()

f_write_03 = open(filename+"_cutoff_n"+".dat", "w")
f_write_03.write("%s\n%d\n%.2f\n%.2f\n%.2f\n%.5f\n%.5f\n%.5f\n%.5f\n%.5f\n%.2f\n%i\n%.2f\n%.2f" % (MAT,len(data),lambda0_in,Tc_in,T_in,t_in,lambdaT_in,f_Eint_0(lambdaT_in),Eint_redttiang,Eint_redttiang/E00,H_appl,par_n,rMax_in,rMax_trian_in))
#MAT=
##Totalvotices=
#lambda0[um]=
#Tc[K]=
#T[K]=
#t=
#lambda(T,Tc)[um]=
#Energy_vortex_line(lambda)[(1e-8)erg/cm=
#Eint_trianglatt(B,lambda)[(1e-8)erg/cm]=
#Eint_trianglatt/E0=
#Happl[G]=
#Parameter=
#rMAX(TF)[um]=
#rMAX(TRIANG_LAT)[um]=
f_write_03.close()

Eint = zeros(len(data))
Eint_E0 = zeros(len(data))
fint_x= zeros(len(data))
fint_y= zeros(len(data))
rij=0
r_ij_in = Vec2d(0,0)
r_i_in = Vec2d(0,0)
r_j_in = Vec2d(0,0)
fint_out = Vec2d(0,0)
fint_mod= Vec2d(0,0)
for i01 in range(len(data)-1):
   for i02 in range(i01+1,len(data)):
	    #rij = sqrt((x_um[i01]-x_um[i02])**2 + (y_um[i01]-y_um[i02])**2)
	    r_i_in = Vec2d(x_um[i01],y_um[i01])
	    r_j_in = Vec2d(x_um[i02],y_um[i02])
	    r_ij_in= r_i_in-r_j_in
	    rij = r_ij_in.get_length()
	    if (rij<0):
			print "Fatal Error:NEGATIVO!!"
	    #print "\n",rij,"\n"
	    if (rij < rMax_in):
	    #if (rij > 0):
			#print "entro"
			Eint[i01] +=f_Eint_K(rij,lambdaT_in)
			Eint[i02] +=f_Eint_K(rij,lambdaT_in)
			fint_out =Vec2d(f_fint_K1(r_i_in,r_j_in,lambdaT_in,E00)[0],f_fint_K1(r_i_in,r_j_in,lambdaT_in,E00)[1])
			fint_x[i01] += fint_out.x
			fint_y[i01] += fint_out.y
			fint_x[i02] -= fint_out.x #fij=-fji
			fint_y[i02] -= fint_out.y
#Eint_E0[i3] = Ein[i3]/f_E12_0(lambdaT_in)
#Normalized energy values
Eint_E0 = Eint/(E00*1.0)
f_write_02 = open(filename+"_Eint_Fint"+".dat", "w")
#for i03 in range(len(data)):
for i03 in range(len(data)):
	fint_mod = Vec2d(fint_x[i03],fint_y[i03])
	f_write_02.write("%.5f %.5f %.6f %.6f %.6f %.6f %.6f\n" % (x_um[i03],y_um[i03],Eint[i03],Eint_E0[i03],fint_x[i03],fint_y[i03],fint_mod.get_length()))
	#print data[indices_exp_out[0][ii],0],data[indices_exp_out[0][ii],1]
f_write_02.close()

#OUTPUT FILE
# x (um) | y (um) | Eint (#OUPUT (1e-8)erg/cm) | Eint/E0 (dimensionless)| f_x | f_y | mod(f)
#forces OUPUT  (1e-4)dyn/cm  = 10^-7 N/m 
