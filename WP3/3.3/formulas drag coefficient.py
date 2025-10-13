import math as m

# Constants
D_fus = 0 #fuselage diameter
L1 = 0 #length first
L2 = 0
L3 = 0
S_wexp = 0 # exposed wing surface area [m^2]
S_hexp = 0 # exposed horizontal tail surface area [m^2]
S_vexp = 0 # exposed vertical tail surface area [m^2]

rho = 0.287 # density [kg/m^3]
V = 200.63 # velocity [m/s]
l = 2.296 # length average root tip chord[m]
mu = 1.8 * 10^-5 # viscosity [Pa*s]
k = 0.052*10**-5 # Surface: smooth molded composite [m]
M = 0.68 # Mach number [-]

xc = 0 #(x/c)m
tc = 0 #thicknes chord airfoil
M = 0.68 #machnumber
sweep = 0 #wingsweep
ffus = (L1+L2+L3)/D_fus
lnac = 0 #length nacele
dnac = 0 #diameter nacele
fnac = lnac/dnac
xct = 0 #(x/c)m for the horizontal tail
sweeptail = 0 #sweep horizontal tail
tct = 0 #thicknes chord airfoil horizontal tail
sweepvtail = 0 #sweep vertical tail
xcv = 0 #(x/c)m vertical tail
tcv = 0 #(t/c) vertical tail

u = 0 # upsweep angle [rad]
Amax = 0 # maximum fuselage cross-sectional area [m^2]
Abase = 0 # fuselage base area [m^2]

SA_LG = 0 # frontal area landing gear [m^2]
d_LG = 0 # height landing gear [m]
w_LG = 0 # width landing gear [m]


#Sw calulation
Swet_wing = 1.07 * 2 * S_wexp
Swet_fus = ((m.pi * D_fus)/4)(1/(3*L1^2)((4*L1^2+D_fus^2/4)^1.5-(D_fus^3/8))-D_fus+4*L2+2*m.sqrt(L3^2+D_fus^2/4))
Swet_hor = 1.05 * 2 * S_hexp
Swet_vert = 1.05 * 2 * S_vexp
Swet_nac = 0 #wetted area nacelle is dependent on team engine


print(Swet_wing)
print(Swet_hor)
print(Swet_vert)
print(Swet_fus)

S_wet = Swet_wing + Swet_fus + Swet_hor + Swet_vert
print(S_wet)


# Cf calculations
Re_actual = rho*V*l/mu
Re_cutoff = 44.62*(l/k)**(1.053)*M**1.16
Re = min(Re_actual, Re_cutoff)

Cf_laminar = 1.328/m.sqrt(Re)
Cf_turbulent = 0.455 / ((m.log10(Re))**2.58*(1+0.144*M**2)**0.65)

Cf_fuselage = 0.1*Cf_laminar + 0.9*Cf_turbulent
Cf_wing = 0.35*Cf_laminar + 0.65*Cf_turbulent
Cf_tail = 0.35*Cf_laminar + 0.65*Cf_turbulent
Cf_nac = 0 #CF nacelle still needed

# FF
FFwing =(1+(0.6/xc)*(tc)+100(tc)^4) (1.34*M^0.18*(m.cos(sweep)^0.28))
FFhort = (1+(0.6/xct)*(tct)+100(tct)^4) (1.34*M^0.18*(m.cos(sweeptail)^0.28))
FFvertc = (1+(0.6/xcv)*(tcv)+100(tcv)^4) (1.34*M^0.18*(m.cos(sweepvtail)^0.28))
FFfus = (1+(60/ffus^3)+(ffus/400))
FFnac = (1+(0.35/fnac))

#IF constants

IFnac = 1.5 #directly under fuselage nacelle
IFwing = 1.2 #low wing config
IFfuselage = 1 #?
IFtail = 1.04 #T-tail accounts for vtail and htail

totalwing = FFwing*Cf_wing*IFwing*Swet_wing
totalfus = FFfus*Cf_fuselage*IFfuselage*Swet_fus
totalvtail = FFvertc*Cf_tail*IFtail*Swet_vert
totalhtail = FFhort*Cf_tail*IFtail*Swet_hor
totalnac = FFnac*Cf_nac*IFnac*Swet_nac

total = totalwing+ totalfus+ totalhtail +totalvtail +totalnac
Sref = 46.919 #wing area


# Miscellaneous drag
Cd_wavedrag = 0 # no wave drag because M < Mcr

Dq_upsweep = 3.83*u**2.5*Amax
Cd_upsweep = Dq_upsweep / Sref # miscellaneous drag due to fuselage upsweep

Dq_base = (0.139 + 0.419*(M - 0.161)**2)*Abase
Cd_base = Dq_base / Sref # fuselage base drag

Ss = d_LG*w_LG # reference area landing gear
C_Ds = 0.05328*m.exp((5.615*SA_LG)/Ss) # for open wheel wells, for closed: 0.04955
Cd_LG = C_Ds*Ss/Sref # landing gear miscellaneous drag

Cdmis = Cd_wavedrag + Cd_upsweep + Cd_base + Cd_LG



#flap constants
Fflap = 0.0074
cfc =  0.25 #flap chord ratio
Sflap = 20.135#area flap
deltaf = 40 #deflection flap

Dflap = Fflap* (cfc)(Sflap/Sref)*(deltaf-10)
# Cd0
Cd0 = 1/Sref * total +Cdmis