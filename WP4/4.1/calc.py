import numpy as np
import scipy as sp

class Calc():
    def __init__(self, file):
        with open(file) as data:
            dat = np.genfromtxt(data, skip_header=21, invalid_raise=False)
            
        ylst = dat[0,:]
        Cllst = dat[3,:]
        Cdlst = dat[5,:]
        Cmlst = dat[7,:]    
        
        self.Cl = sp.interpolate.interpld(ylst, Cllst, kind='linear', fill_value='extrapolate')
    
    # Load Distribution as function of x
    def distrib(self):
        distrib = None
        return distrib
    
    # Normal force as function of x
    def axial(self, x, L, distrib):
        N = sp.integrate.quad(distrib, x, L) 
        return N

    # Shear force as function of x
    def shear(self, x, L, distrib):
        V = sp.integrate.quad(distrib, x, L) 
        return V