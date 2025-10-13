import math as m

# Constants
D_fus = 0 #fuselage diameter
L1 = 0 #length first
L2 = 0
L3 = 0
S_wexp = 0
S_hexp = 0
S_vexp = 0

rho = 0
V = 0
l = 0
mu = 0
k = 0
M = 0

xc = 0 #(x/c)m
tc = 0 #thicknes chord airfoil
M = 0 #machnumber
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


# FF
FFwing =(1+(0.6/xc)*(tc)+100(tc)^4) (1.34*M^0.18*(m.cos(sweep)^0.28))
FFhort = (1+(0.6/xct)*(tct)+100(tct)^4) (1.34*M^0.18*(m.cos(sweeptail)^0.28))
FFvertc = (1+(0.6/xcv)*(tcv)+100(tcv)^4) (1.34*M^0.18*(m.cos(sweepvtail)^0.28))
FFfus = (1+(60/ffus^3)+(ffus/400))
FFnac = (1+(0.35/fnac))