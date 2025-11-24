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
        
    # Normal force as function of x
    def axial(self, x, L, loading, pointLoads):
        N = sp.integrate.quad(loading, x, L) 
        
        for i in range(len(pointLoads)):
            N += pointLoads[i,1] * (1-np.heaviside(x-pointLoads[i,0], 1))
        
        return N(x)

    # Shear force as function of x
    def shear(self, x, L, loading, pointLoads):
        V = sp.integrate.quad(loading, x, L) 
        return V

    # Moment force as function of x
    def moment(self, x, L, distrib, M_p, u):
        u_p = u / L
        M = -(sp.integrate.quad(distrib, x, L) + M_p * (1 - u_p))
        return M

    def torque(self, ):
        T = sp.integrate.quad()
        return T

