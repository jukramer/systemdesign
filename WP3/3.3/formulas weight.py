from parameters import *
import math

# Wing weight estimation
W_W = 0.0017*W_MZF*(b/math.cos(Lambda_12))**0.75*(1 + (6.3*math.cos(Lambda_12)/b)**(1/2))*n_ult**0.55*(b*S/(t_r*W_MZF*math.cos(Lambda_12)))**0.30

# Empennage weight estimation
W_h = k_h*S_h*(3.81*(S_h**0.2*V_D)/((1000*math.cos(Lambda_12_h))**(1/2))-0.287)
W_v = k_v*S_v*(3.81*(S_v**0.2*V_D)/((1000*math.cos(Lambda_12_v))**(1/2))-0.287)

W_emp = W_h + W_v

# Fuselage weight estimation



#Flight control system
W_fc = K_fc * (W_TO)^(2/3)

#Hydraulic systems
W_hy = 0.008 * W_TO


#electrical systems
W_ELS = 10.8(V_pax)^0.7*(1-0.018(V_pax)^0.35)

#instrumentation, avionics and electronics
W_IAE = 0.575(W_E)^0.556 * R^0.25

#airconditioning pressurization, ani- and- deicing systems
W_api = 6.75(l_pax)^1.28

#oxygen system
W_ox = 40 + 2.4* N_pax

#auxilary power unit
W_apu = 0.08 * W_TO


