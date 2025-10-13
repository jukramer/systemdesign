import math

rho = 0
V = 0
l = 0
mu = 0
k = 0
M = 0


# Cf calculations
Re_actual = rho*V*l/mu
Re_cutoff = 44.62*(l/k)**(1.053)*M**1.16
Re = min(Re_actual, Re_cutoff)

Cf_laminar = 1.328/math.sqrt(Re)
Cf_turbulent = 0.455 / ((math.log10(Re))**2.58*(1+0.144*M**2)**0.65)

Cf_fuselage = 0.1*Cf_laminar + 0.9*Cf_turbulent
Cf_wing = 0.35*Cf_laminar + 0.65*Cf_turbulent
Cf_tail = 0.35*Cf_laminar + 0.65*Cf_turbulent

# Misc