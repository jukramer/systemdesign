import numpy as np
import matplotlib.pyplot as plt
from parameters import *

class Calc():
    def __init__(self):
        pass
    
    def v_stall(self, W, S, rho, C_Lmax):
        WS = W/S
        return np.sqrt(WS*2/rho/C_Lmax)
    
    def v_max(self, T, W, S, rho, C_D0, A, e):
        T_alt = T*(rho/rho_SL)
        TWmax = T_alt/W
        WS = W/S
        
        return np.sqrt((TWmax*WS+WS*np.sqrt(TWmax**2-4*C_D0/(np.pi*e*A)))/(rho*C_D0))
    
    def h_to_rho(self, h_vals):
        pass
    