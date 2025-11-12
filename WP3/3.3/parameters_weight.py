# Weight
m_to_ft = 3.280839895           # meter to feet conversion
g = 9.80665 * m_to_ft           # gravitational acceleration [ft/s^2]
lbs_to_N = 4.448                # pound to newton conversion
kg_to_lbs = 2.20462262          # kg to lbs

W_TO = 9907.7*kg_to_lbs         # maximum take off weight [lbs]
W_F = 2960*kg_to_lbs            # fuel weight [lbs]
W_MZF = W_TO - W_F              # maximum zero fuel weight [lbs]
b = 20.435*m_to_ft              # span [ft]
Lambda_12 = 0                   # wing semi-chord sweep angle [rad]
S = 46.919*m_to_ft**2           # wing area [ft^2]
n_ult = 4.275                   # ultimate load factor [-]
t_r = 0                         # maximum root thickness [ft]

k_h = 1.0                       # for fixed incidence stabilizers [-]
S_h = 0                         # horizontal tail area [ft^2]
V_D = 0                         # design dive speed [KEAS]
Lambda_12_h = 0                 # semi-chord sweep angle horizontal tail [rad]

S_v = 11.37711*m_to_ft**2                         # vertical tail area [ft^2]
z_h = 0                         # distance from vert.tail root to where ht is mounted on the v.t. [ft]
b_v = 0                         # vertical tail span [ft]
k_v = 1 + 0.15*(S_h*z_h/(S_v*b_v))         # for fuselage mounted vertical tails
Lambda_12_v = 0                 # semi-chord sweep angle vertical tail [rad]

k_f = 1.08          # for a pressurized fuselage
V_D = 0             # design dive speed [KEAS]
l_h = 0             # distance from wing root c/4 to hor. tail root c/4 [ft]
w_f = 0             # maximum fuselage width [ft]
h_f = 0             # maximum fuselage height [ft]
S_fgs = 0           # fuselage gross shell area [ft]

T_TO = 0            # total required take-off thrust

K_gr = 1.0          # for low wing airplanes [-]
A_g_main = 33.0     # constant for land gear main [-]
B_g_main = 0.04     # constant for land gear main [-]
C_g_main = 0.021    # constant for land gear main [-]
D_g_main = 0.0      # constant for land gear main [-]
A_g_nose = 12.0     # constant for land gear nose [-]
B_g_nose = 0.06     # constant for land gear nose [-]
C_g_nose = 0.0      # constant for land gear nose [-]
D_g_nose = 0.0      # constant for land gear nose [-]
K_fc = 0.64         # constant for airplanes with powerd flight controls [-]
V_pax = 0           # passenger cabin volime in [ft^3]
R = 0               # range in nautical miles
W_E = 0             # empty weight in [lbs]
l_pax = 0           # length of the passenger cabin [ft]
N_pax = 6           # number of passengers



c = 10 #wing mean geometric chord [ft]
CL_alpha = 0.5 #CL alpha
h_cr = 41000 #cruise altitude [ft]
Gw = 90 #flight design gross weight [lbs]
W_cr = 100 #cruise weight [lbs]
CL_max = 0 