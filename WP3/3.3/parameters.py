# Drag coefficient
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



# Weight
g = 9.80665                     # gravitational acceleration [m/s^2]
lbs_to_N = 4.448                # pound to newton conversion
m_to_ft = 3.280839895           # meter to feet conversion
W_TO = 9907.7*g                 # maximum take off weight [N]
W_F = 2960*g                    # fuel weight [N]
W_MZF = (W_TO - W_F)/lbs_to_N   # maximum zero fuel weight
b = 20.435*m_to_ft              # span [ft]
Lambda_12 = 0                   # wing semi-chord sweep angle [rad]
S = 46.919*m_to_ft**2           # wing area [ft^2]
n_ult = 4.275                   # ultimate load factor [-]
t_r = 0                         # maximum root thickness [ft]

k_h = 1.0                       # for fixed incidence stabilizers [-]
S_h = 0                         # horizontal tail area [ft^2]
V_D = 0                         # design dive speed [KEAS]
Lambda_12_h = 0                 # semi-chord sweep angle horizontal tail [rad]

S_v = 0                         # vertical tail area [ft^2]
z_h = 0                         # distance from vert.tail root to where ht is mounted on the v.t. [ft]
b_v = 0                         # vertical tail span [ft]
k_v = 1 + 0.15*(S_h*z_h/(S_v*b_v))         # for fuselage mounted horizontal tails
Lambda_12_v = 0                 # semi-chord sweep angle horizontal tail [rad]

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
W_E = 0             # empty weight in[lbs]
l_pax = 0           # length of the passenger cabin [ft]
N_pax = 6           # number of passengers

