from math import pi

# Drag coefficient
D_fus = 2.72 #fuselage diameter
L1 = 5.445 #length first
L2 = 3.956  
L3 = 8.712
S_wexp = 38 # exposed wing surface area [m^2]
S_hexp = 0 # exposed horizontal tail surface area [m^2]
S_vexp = 0 # exposed vertical tail surface area [m^2]

rho = 0.287 # density [kg/m^3]
V = 200.63 # velocity [m/s]
l = 2.296 # length average root tip chord[m]
mu = 1.8 * 10**-5 # viscosity [Pa*s]
k = 0.052*10**-5 # Surface: smooth molded composite [m]
M = 0.68 # Mach number [-]

xc = 0 #(x/c)m
tc = 0.12 #thicknes chord airfoil
M = 0.68 #machnumber
sweep = 10.56*(pi/180) #wingsweep
ffus = (L1+L2+L3)/D_fus
lnac = 2.502 #length nacele
dnac = 1.35 #diameter nacele
fnac = lnac/dnac
xct = 0 #(x/c)m for the horizontal tail
sweeptail = 25*(pi/180) #sweep horizontal tail
tct = 0 #thicknes chord airfoil horizontal tail
sweepvtail = 35*(pi/180) #sweep vertical tail
xcv = 0 #(x/c)m vertical tail
tcv = 0 #(t/c) vertical tail

u = 5*(pi/180) # upsweep angle [rad]
Amax = 2.722 # maximum fuselage cross-sectional area [m^2]
#Abase = 0 # fuselage base area [m^2]

SA_LG = 0 # frontal area landing gear [m^2]
d_LG = 0 # height landing gear [m]
w_LG = 0 # width landing gear [m]



