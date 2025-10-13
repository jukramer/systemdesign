import numpy as np
import matplotlib.pyplot as plt
from parameters import *

class Calc():
    def __init__(self):
        pass
    
    def v_stall(self, WS, rho, C_Lmax):
        return np.sqrt(WS*2/rho/C_Lmax)
    
    def v_max(self, T, WS, rho, C_D0, A, e):
        T_alt = T*(rho/rho_SL)
        W = WS*S
        
        TWmax = T_alt/W
        
        return np.sqrt((TWmax*WS+WS*np.sqrt(TWmax**2-4*C_D0/(np.pi*e*A)))/(rho*C_D0))
    
    def h_to_rho(self, h_vals):
        return self.ISA(h_vals)[2]
    
    def ISA(self, h):
        g = 9.80665
        R = 287
        T0 = 288.15
        p0 = 101325
        rho = 1.225
        h0 = 0

        layers = [
            (11000, -0.0065),
            (20000, 0),
            (32000, 0.001),
            (47000, 0.0028),
            (51000, 0),
            (71000, -0.0028),
            (86000, -0.002)]

        for hLimit, a in layers:
            if h <= h0:
                break

            h1 = min(h, hLimit)
            delta_h = h1-h0

            if a == 0:
                T1 = T0
                p1 = p0 * np.exp((-g*delta_h/(R*T0)))
            else:
                T1 = T0 + a * delta_h
                p1 = p0 * (T1/T0)**(-g/(a*R))

            T0 = T1
            p0 = p1
            h0 = h1
            
        density = p0/(R*T0)
        
        return T0, p0, density 
    
    def plotFE(self, h_min, h_max, h_step):
        h_vals = np.arange(h_min, h_max, h_step)
        T_vals = []
        rho_vals = []
        for h in h_vals:
            T_vals.append(self.ISA(h)[0])
            rho_vals.append(self.ISA(h)[2])
            
        T_vals = np.array(T_vals)
        rho_vals = np.array(rho_vals)
        
        v_stall_vals = self.v_stall(WS_CR, rho_vals, CL_MAX)
        v_max_vals = self.v_max(TMAX_SL, WS_CR, rho_vals, C_D0, ASPECT_RATIO, e)
        
        M_stall_vals = v_stall_vals/np.sqrt(GAMMA*R_AIR*T_vals)
        M_max_vals = v_max_vals/np.sqrt(GAMMA*R_AIR*T_vals)
        
        plt.plot(M_stall_vals, h_vals/0.3048)
        plt.plot(M_max_vals, h_vals/0.3048)
        plt.show()