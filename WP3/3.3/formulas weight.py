from parameters import *
import math

# Wing weight estimation
W_W = 0.0017*W_MZF*(b/math.cos(Lambda_12))**0.75*(1 + (6.3*math.cos(Lambda_12)/b)**(1/2))*n_ult**0.55*(b*S/(t_r*W_MZF*math.cos(Lambda_12)))**0.30

# Empennage weight estimation
W_h = k_h*S_h*(3.81*(S_h**0.2*V_D)/((1000*math.cos(Lambda_12_h))**(1/2))-0.287)
W_v = k_v*S_v*(3.81*(S_v**0.2*V_D)/((1000*math.cos(Lambda_12_v))**(1/2))-0.287)

W_emp = W_h + W_v

# Fuselage weight estimation
W_f = 0.021*k_f*(V_D*l_h/(w_f + h_f))^(1/2)*S_fgs**1.2

# Nacelle weight estimation (high bypass ratio turbofan)
W_n = 0.065*T_TO

# Landing gear weight estimation
# main landing gear
W_g_main = K_gr*(A_g_main + B_g_main*(W_TO)**(3/4) + C_g_main*W_TO + D_g_main*(W_TO)**(3/2))

# nose landing gear
W_g_nose = K_gr*(A_g_nose + B_g_nose*(W_TO)**(3/4) + C_g_nose*W_TO + D_g_nose*(W_TO)**(3/2))

W_g = W_g_main + W_g_nose

# Total structure weight
W_struc = W_W + W_emp + W_f + W_n + W_g     # [lbs]
W_struc = lbs_to_N*W_struc                     # [N]