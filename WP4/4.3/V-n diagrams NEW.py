import math
import numpy as np
import matplotlib.pyplot as plt

#print(V_s0_approach)
print(V_s0_landing)
#print(V_s0_takeoff)


def case1():
    # Weights / parameters
    W_payload = 1010 # kg
    W_oem = 4876.78 # kg
    W_mtom = 7142.310 # kg
    W_mtom_lbs = 0.454 * W_mtom # lbs
    S = 33.609 # m^2
    C_L_max_cr = 1.663
    C_L_max_to = 2.1
    C_L_max_a = 2.1
    C_L_max_l = 2.1
    
    # Flight conditions
    M_cr = 0.68
    Altitude = 41000 # ft
    Altitude_metric = Altitude / 3.281 # m
    rho_cr = 0.287407 # kg/m3 or m/s conversion factor
    g = 9.80665 # m/s^2
    rho = 1.225 # kg/m^3
    pressure_scaling = math.sqrt(rho/rho_cr)
    
    # Stall speeds
    V_s0_takeoff = np.sqrt(2*g*W_mtom/(rho*S*C_L_max_to))   # flaps down takeoff
    V_s0_approach = np.sqrt(2*g*W_mtom/(rho*S*C_L_max_a))   # flaps down approach
    V_s0_landing = np.sqrt(2*g*W_mtom/(rho*S*C_L_max_l))   # flaps down landing
    V_s1 = np.sqrt(2*g*W_oem/(rho_cr*S*C_L_max_cr))  * pressure_scaling # flaps up
    
    # Limit load factors
    n_max_flaps = 2
    n_max = 3.192
    n_max = max(2.5, min(n_max, 3.8))
    n_min = -1
    
    # V-speeds
    V_c = (M_cr * rho_cr) * pressure_scaling
    V_d = (V_c / 0.8) * pressure_scaling
    V_a = (V_s1 * math.sqrt(n_max)) * pressure_scaling
    
    # Flap speeds
    V_f_takeoff = 1.6 * V_s1
    V_f_approach = 1.8 * V_s1
    V_f_landing = 1.8 * V_s0_landing
    
    return {
        'V_s0_takeoff': V_s0_takeoff, 'V_s0_approach': V_s0_approach, 'V_s0_landing': V_s0_landing, 'V_s1': V_s1, 'V_a': V_a, 'V_c': V_c, 'V_d': V_d,
        'V_f_takeoff': V_f_takeoff, 'V_f_approach': V_f_approach, 'V_f_landing': V_f_landing,
        'n_max': n_max, 'n_min': n_min, 'n_max_flaps': n_max_flaps
    }

########## GRAPHING #############

case_data = case1()
V_d = case_data['V_d']
V_s0_takeoff = case_data['V_s0_takeoff']
V_s0_approach = case_data['V_s0_approach']
V_s0_landing = case_data['V_s0_landing']
V_s1 = case_data['V_s1']
V_a = case_data['V_a']
V_c = case_data['V_c']
n_max = case_data['n_max']
n_min = case_data['n_min']
V_f_takeoff = case_data['V_f_takeoff']
V_f_approach = case_data['V_f_approach']
V_f_landing = case_data['V_f_landing']


n_max_flaps = 2

def line_n_pos(V):
    if V <= V_a:
        n_pos = (V / V_s1)**2
    else:
        n_pos = n_max
    return n_pos

def line_n_flaps(V):
    if V <= V_a:
        n_pos = (V / V_s1)**2
    else:
        n_pos = n_max
        
    if V <= V_f_takeoff:
        n_flaps_calc = (V / V_s0_takeoff)**2
        if n_flaps_calc <= n_max_flaps:
            n_flaps = n_flaps_calc
        else:
            n_flaps = n_max_flaps
    else:
        n_flaps = float('nan')
        
    if n_flaps <= n_pos:
        n_flaps = float('nan')
    
    return n_flaps

def line_n_neg(V):
    n_neg_calc = -(V / V_s1)**2
    if n_neg_calc >= -1:
        n_neg = n_neg_calc
    elif V <= V_c:
        n_neg = -1

    else:
        n_neg = -1 + (V - V_c) / (V_d - V_c)
    return n_neg

V_values = np.linspace(0, V_d, 200)

n_pos_values = []
n_neg_values = [] 
n_flaps_values = []

for V_val in V_values:
    n_pos = line_n_pos(V_val)
    n_neg_val = line_n_neg(V_val)
    n_flaps_val = line_n_flaps(V_val)
    n_pos_values.append(n_pos)
    n_neg_values.append(n_neg_val)
    n_flaps_values.append(n_flaps_val)

########## PLOTTING #############
'''
plt.figure(figsize=(10, 6))
plt.xlim(0, V_d; + 10)  # Set X-axis limits
plt.ylim(n_min - 0.5, n_max + 0.5)  # Set Y-axis limits
plt.plot([0, V_d + 10], [0, 0], 'k-', linewidth=1.5)
plt.plot(V_values, n_pos_values, 'b-', linewidth=2)
plt.plot(V_values, n_neg_values, 'b-', linewidth=2) 
plt.plot(V_values, n_flaps_values, 'b-', linewidth=2)
plt.plot([V_d, V_d], [0, n_max], 'b-', linewidth=2)
plt.plot([V_a, V_a], [0, n_max], 'k--', linewidth=1)
plt.plot([V_c, V_c], [n_min, n_max], 'k--', linewidth=1)
plt.plot([0, V_d], [1, 1], 'k--', linewidth=1)
plt.plot([0, V_d], [-1, -1], 'k--', linewidth=1)
plt.plot([V_s1, V_s1], [0, 1], 'k--', linewidth=1)
plt.title('V_EAS vs n Diagram')
plt.xlabel('V_EAS')
plt.ylabel('n')
#plt.legend()
plt.grid(True)
plt.axis
plt.show()'''


fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1 
axes[0].set_xlim(0, V_d + 10)
axes[0].set_ylim(n_min - 0.5, n_max + 0.5)
axes[0].plot([0, V_d + 10], [0, 0], 'k-', linewidth=1.5)
axes[0].plot(V_values, n_pos_values, 'b-', linewidth=2)
axes[0].plot(V_values, n_neg_values, 'b-', linewidth=2) 
axes[0].plot(V_values, n_flaps_values, 'b-', linewidth=2)
axes[0].plot([V_d, V_d], [0, n_max], 'b-', linewidth=2)
axes[0].plot([V_a, V_a], [0, n_max], 'k--', linewidth=1)
axes[0].plot([V_c, V_c], [n_min, n_max], 'k--', linewidth=1)
axes[0].plot([0, V_d], [1, 1], 'k--', linewidth=1)
axes[0].plot([0, V_d], [-1, -1], 'k--', linewidth=1)
axes[0].plot([V_s1, V_s1], [0, 1], 'k--', linewidth=1)
axes[0].set_title('V_EAS vs n Diagram, Flaps for Take-off, MTOM')
axes[0].set_xlabel('V_EAS')
axes[0].set_ylabel('n')
#axes[0].legend()
axes[0].grid(True)

# Plot 2
axes[1].set_xlim(0, V_d + 10)
axes[1].set_ylim(n_min - 0.5, n_max + 0.5)
axes[1].plot([0, V_d + 10], [0, 0], 'k-', linewidth=1.5)
axes[1].plot(V_values, n_pos_values, 'b-', linewidth=2)
axes[1].plot(V_values, n_neg_values, 'b-', linewidth=2) 
axes[1].plot(V_values, n_flaps_values, 'b-', linewidth=2)
axes[1].plot([V_d, V_d], [0, n_max], 'b-', linewidth=2)
axes[1].plot([V_a, V_a], [0, n_max], 'k--', linewidth=1)
axes[1].plot([V_c, V_c], [n_min, n_max], 'k--', linewidth=1)
axes[1].plot([0, V_d], [1, 1], 'k--', linewidth=1)
axes[1].plot([0, V_d], [-1, -1], 'k--', linewidth=1)
axes[1].plot([V_s1, V_s1], [0, 1], 'k--', linewidth=1)
axes[1].set_title('V_EAS vs n Diagram, Flaps for Approach, MTOM')
axes[1].set_xlabel('V_EAS')
axes[1].set_ylabel('n')
#axes[1].legend()
axes[1].grid(True)

# Plot 3
axes[2].set_xlim(0, V_d + 10)
axes[2].set_ylim(n_min - 0.5, n_max + 0.5)
axes[2].plot([0, V_d + 10], [0, 0], 'k-', linewidth=1.5)
axes[2].plot(V_values, n_pos_values, 'b-', linewidth=2)
axes[2].plot(V_values, n_neg_values, 'b-', linewidth=2) 
axes[2].plot(V_values, n_flaps_values, 'b-', linewidth=2)
axes[2].plot([V_d, V_d], [0, n_max], 'b-', linewidth=2)
axes[2].plot([V_a, V_a], [0, n_max], 'k--', linewidth=1)
axes[2].plot([V_c, V_c], [n_min, n_max], 'k--', linewidth=1)
axes[2].plot([0, V_d], [1, 1], 'k--', linewidth=1)
axes[2].plot([0, V_d], [-1, -1], 'k--', linewidth=1)
axes[2].plot([V_s1, V_s1], [0, 1], 'k--', linewidth=1)
axes[2].set_title('V_EAS vs n Diagram, Flaps for Landing, MTOM')
axes[2].set_xlabel('V_EAS')
axes[2].set_ylabel('n')
#axes[2].legend()
axes[2].grid(True)

# Common elements
for ax in axes:
    ax.plot([0, V_d + 10], [0, 0], 'k-', linewidth=1.5)
    ax.set_xlabel('V_EAS')
    ax.set_ylabel('n')

plt.tight_layout()
plt.show()
