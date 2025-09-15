import numpy as np
from parameters import *

class Calc():
    def __init__(self, INT_INTRVLS):
        self.INT_INTRVLS = INT_INTRVLS
    
    def c(self, y):
        return c_ROOT - (c_ROOT - c_TIP)*(2*y/(b))
    
    def mac(self):
        y_VALS = np.arange(0, b/2, self.INT_INTRVLS)
        return 2/S*np.trapz(self.c(y_VALS), y_VALS)
    
    def pos_mac(self):
        y = (b/6)*((1+2*TAPER_RATIO)/(1+TAPER_RATIO))
        x = y*np.tan(LAMBDA_LE)
        
        return x, y
    
    def 
    