# REQUIREMENTS
M_MAX_PL= 1010 # [kg]
MACH_CRUISE = 0.68 # @ 41000ft

TAKEOFF_DIST = 1250 # [m]
LANDING_DIST = 620 # [m]

RANGE_DESIGN = 4482e3 # [m] @ 1010kg
RANGE_HARMONIC = 6.1e6 # [m] @ 1010kg
RANGE_FERRY = 7e6 # [m] 

# AIRCRAFT PROPERTIES - FOR NOW ASSUME FAR/25
M_OEW = 10000 # [kg]
M_MTOW = 12000 # [kg]

CL_MAX_CR = 1.6 # Cruise
CL_MAX_TO = 1.9 # Takeoff
CL_MAX_L = 2.1 # Landing
AR = 5 # Aspect ratio
CD_0 = 0.2 # CD0

C_LFL = 0.45 # FAR 25 0.6 for 23
V_APP = 60 # [m/s] estimate

B = 5 # Bypass ratio, assumed
THETA_BREAK = 1.07 # assumed

# ENVIRONMENTAL CONDITIONS
RHO_SL = 1.225 # [kg/m^3]
GAMMA = 1.4 # specific heat ratio in air

n_v = 0.85

# ENVIRONMENTAL CONDITIONS AT CRUISE ALTITUDE
alt_cr = 41000 #ft
    = 12496.8 #m
Temp_cr = 216.65 #K
p_cr = 17868.132 #Pa
rho_cr = 0.287368 #kg/m^3
T_SL = 288.15 # [K]
P_SL = 101325 # [Pa]
T_CR = 200
P_CR = 50000


