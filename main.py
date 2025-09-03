import matplotlib.pyplot as plt
import numpy as np
from parameters import *

L = 1/2 * rho * V ** 2 * C_l * S #lift equation
V_app = 1.23*V_so # relation between approach and stall speed in CS-25 plane
N_f = m_f/(m_f+m_t) # gravimetric efficiency
e_eff = η_f * rho_e/rho # effective specific energy for fuels
rho_e_eff = n_v * rho_e # effective energy density for fuels
