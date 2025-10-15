from parameters import *
import math

# Wing weight estimation
W_W = 0.0017*W_MZF*(b/math.cos(Lambda_12))**0.75*(1 + (6.3*math.cos(Lambda_12)/b)**(1/2))*n_ult**0.55*(b*S/(t_r*W_MZF*math.cos(Lambda_12)))**0.30

# Empennage weight estimation
W_h = k_h*S_h*(3.81*(S_h**0.2*V_D)/((1000*math.cos(Lambda_12_h))**(1/2))-0.287)
W_v = k_v*S_v*(3.81*(S_v**0.2*V_D)/((1000*math.cos(Lambda_12_v))**(1/2))-0.287)

W_emp = W_h + W_v

# Fuselage weight estimation

