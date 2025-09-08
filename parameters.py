# REQUIREMENTS
M_MAX_PL= 1010 # [kg]
MACH_CRUISE = 0.68 # @ 41000ft

TAKEOFF_DIST = 1250 # [m]
LANDING_DIST = 620 # [m]

RANGE_DESIGN = 4482e3 # [m] @ 1010kg
RANGE_HARMONIC = 6.1e6 # [m] @ 1010kg
RANGE_FERRY = 7e6 # [m] 

# AIRCRAFT PROPERTIES - FOR NOW ASSUME FAR/25
M_OEW = 7000 # [kg] DUMMY
M_MTOW = 11000 # [kg] DUMMY

CL_MAX_CR = 1.6 # Cruise
CL_MAX_TO = 1.9 # Takeoff
CL_MAX_L = 2.1 # Landing
AR = 8.5 # Aspect ratio
CD_0 = 0.0228 # CD0

C_LFL = 0.45 # FAR 25: 0.45, FAR 23: 0.6
V_APP = 60 # [m/s] estimate

B = 5 # Bypass ratio, assumed
THETA_BREAK = 1.07 # assumed
N_E = 2 # number of engines
KT = 0.85 # from adsee reader

H2 = 11 # from adsee reader, for normal cs25 ac

BETA_TO = 1 # Takeoff mass fraction
BETA_CR = 1 # Cruise mass fraction
BETA_L = 0.85 # Landing mass fraction

# ENVIRONMENTAL CONDITIONS
GAMMA = 1.4 # specific heat ratio in air
R = 287 # gas constant in air

alt_cr = 41000 #ft = 12496.8m
T_SL = 288.15 # [K]
P_SL = 101325 # [Pa]
RHO_SL = 1.225 # [kg/m^3]

T_CR = 216.65 # [K]
P_CR = 17868.132 # [Pa]
RHO_CR = 0.287368 # [kg/m^3]

e = 0.812  #OSWALD EFFICIENCY FACTOR

#CLIMB GRADIENTS
T_CG = 288.15 # [K]
P_CG = 101325 # [Pa]
RHO_SL = 1.225 # [kg/m^3]

#CLIMB RATE
c = 18 # [m/s]
T_CR = 216.65 # [K]
P_CR = 17868.132 # [Pa]
RHO_CR = 0.287368 # [kg/m**3]


