import numpy as np
from parameters import *

class Calc():
    def __init__(self, INT_INTRVLS):
        self.INT_INTRVLS = INT_INTRVLS
    
    def lambda_c4_to_le(self):
        return np.rad2deg(np.arctan(np.tan(np.deg2rad(LAMBDA_C4))-C_ROOT/(2*b)*(TAPER_RATIO-1))) 
    
    def calc_c(self):
        C_ROOT = 2*S/(b*1+TAPER_RATIO)
        print(C_ROOT, TAPER_RATIO*C_ROOT)
    
    def c(self, y):
        return C_ROOT - (C_ROOT - C_TIP)*(2*y/(b))
    
    def mac(self):
        y_VALS = np.arange(0, b/2, self.INT_INTRVLS)
        return 2/S*np.trapz(self.c(y_VALS), y_VALS)
    
    def pos_mac(self):
        y = (b/6)*((1+2*TAPER_RATIO)/(1+TAPER_RATIO))
        x = y*np.tan(LAMBDA_LE)
        
        return x, y
    
    
    def CL_Max(self, airfoil): # airfoil should be 1, 2, 3, 4, or 5, see AIRFOILS in parameters.py
        DELTA_Y = AIRFOILS[airfoil+1][1]
        delta_Y = TC_MAX*DELTA_Y
        
        # Use graph to find CL_MAX
        
    def alpha_stall(self):
        
    
    
    
    
    
    
    
    
    