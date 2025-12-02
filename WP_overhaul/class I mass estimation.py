import matplotlib.pyplot as plt
import numpy as np
import colorsys as cs
from parameters import *




CD_0 = C_f * S_wet_to_S_ratio

e = 1/(np.pi * aspect_ratio * psi + 1/phi)

lift_to_drag_ratio  = 1/2 * np.sqrt(np.pi * aspect_ratio * e / CD_0)

#unfinished
V_cr = M_cr * np.sqrt(gamma * R * T_cr)  # Cruise speed [m/s]

range_lost = 1/0.7 * lift_to_drag_ratio * (h_cr + V_cr**2/(2*g)) # Lost range [m]

range_eq = (range_design + range_lost) * (1 + f_margin) + 1.2*range_diversion + t_l * V_cr

TSFC = 22 * bypass_ratio**(-0.19) * 1e-6 # Thrust specific fuel consumption [1/s]

jet_efficiency = V_cr / (TSFC * e_fuel) # Jet efficiency [-]

mass_fraction_fuel_MTOM =  1 - np.exp(-range_eq/(jet_efficiency * (e_fuel/g) * lift_to_drag_ratio))

m_MTOM = m_payload / (1 - mass_fraction_OEM_MTOM - mass_fraction_fuel_MTOM)

m_OEM = m_MTOM * mass_fraction_OEM_MTOM



print("CD_0: ", CD_0)
print("e: ", e)
print("Lift to drag ratio: ", lift_to_drag_ratio)
print("Cruise speed: ", V_cr)
print("Range lost: ", range_lost)
print("Range equivalent: ", range_eq)
print("TSFC: ", TSFC)
print("Jet efficiency: ", jet_efficiency)
print("Mass fraction fuel MTOM: ", mass_fraction_fuel_MTOM)
print("MTOM: ", m_MTOM)
print("OEM mass: ", m_OEM)
