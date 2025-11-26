from parameters import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

# CONSTANTS
NULL_ARRAY_2 = np.zeros((2,1)) # 0-load array (for 2-row point loads)
NULL_ARRAY_3 = np.zeros((3,1)) # 0-load araray (for 3-row point loads)


class Calc():
    def __init__(self, file0, file10):
        with open(file0) as data0:
            dat0 = np.genfromtxt(data0, skip_header=21, invalid_raise=False)

        with open(file10) as data10:
            dat10 = np.genfromtxt(data10, skip_header=21, invalid_raise=False)

        ylst0 = dat0[:,0]
        Cllst0 = dat0[:,3]
        Cdlst0 = dat0[:,5]
        Cmlst0 = dat0[:,7]    
        clst0 = dat0[:,1]

        ylst10 = dat10[:,0]
        Cllst10 = dat10[:,3]
        Cdlst10 = dat10[:,5]
        Cmlst10 = dat10[:,7]    
        clst10 = dat10[:,1]

        # Interpolate datapoints from XFLR5 data to obtain python functions
        self.Cl0 = sp.interpolate.interp1d(ylst0, Cllst0, kind='cubic', fill_value='extrapolate')
        self.Cd0 = sp.interpolate.interp1d(ylst0, Cdlst0, kind='cubic', fill_value='extrapolate')
        self.Cm0 = sp.interpolate.interp1d(ylst0, Cmlst0, kind='cubic', fill_value='extrapolate')

        self.Cl10 = sp.interpolate.interp1d(ylst10, Cllst10, kind='cubic', fill_value='extrapolate')
        self.Cd10 = sp.interpolate.interp1d(ylst10, Cdlst10, kind='cubic', fill_value='extrapolate')
        self.Cm10 = sp.interpolate.interp1d(ylst10, Cmlst10, kind='cubic', fill_value='extrapolate')

        self.lifttest()
    
    def lifttest(self):
        L =W_MTOW*n_ult
        C_L = L/(q*S)
        print(C_L)

    


    def chord(self, y):
        c = C_ROOT+(C_TIP-C_ROOT)*2*y/b
        return c

    def liftUnitSpan(self, y):
        L = self.Cl(y)*q*self.chord(y)
        return L
    
    def dragUnitSpan(self, y):
        D = self.Cd(y)*q*self.chord(y)
        return D
    
    def momentUnitSpan(self, y):
        M = self.Cm(y)*q*self.chord(y)**2
        return M
    
    def inertialLoading(self, y, massWing, n=1):
        weightDens = n*g*massWing/S
        
        w = weightDens*self.c(y)
        
        return w   
        
    # Normal force as function of x. pointLoads must have cols (position, load) (shape 2xn).
    # loading must be a python function.
    def axial(self, x, L, loading, pointLoads):
        # These asserts ensure that the point load arrays have the correct dimensions
        assert pointLoads.shape[0] == 2
        N = sp.integrate.quad(loading, x, L)[0]
        
        for i in range(pointLoads.shape[1]):
            N += pointLoads[1,i] * (1-np.heaviside(x-pointLoads[0,i], 1))
        
        return N

    # Shear force as function of x. pointLoads must have cols (position, load) (shape 2xn)
    # loading must be a python function.
    def shear(self, x, L, loading, pointLoads):
        # These asserts ensure that the point load arrays have the correct dimensions
        assert pointLoads.shape[0] == 2
        
        V = sp.integrate.quad(loading, x, L)[0]

        for i in range(pointLoads.shape[1]):
            V += pointLoads[1,i] * (1-np.heaviside(x-pointLoads[0,i], 1))
        
        return V

    # Moment force as function of x. pointLoads must have cols (position, load) (shape 2xn)
    # loading must be a python function.
    def moment(self, x, L, loading, pointLoads, args):
        # These asserts ensure that the point load arrays have the correct dimensions
        assert pointLoads.shape[0] == 2
        M = sp.integrate.quad(loading, x, L, args)[0]
        
        for i in range(pointLoads.shape[1]):
            M += pointLoads[1,i] * (1-np.heaviside(x-pointLoads[0,i], 1))
            
        return -M

    # Torsion as function of x. forcePointLoads must have cols (position, load, dist) (shape 3xn);
    # torquePointLoads must have cols (position, load) (shape 2xn)
    # forceloading, torsionLoading, and loadingDist must be a python function.
    def torsion(self, x, L, forceLoading, loadingDist, torsionLoading, forcePointLoads, torquePointLoads):
        # These asserts ensure that the point load arrays have the correct dimensions
        assert forcePointLoads.shape[0] == 3
        assert torquePointLoads.shape[0] == 2
       
        T = sp.integrate.quad(lambda x : forceLoading(x)*loadingDist(x) + torsionLoading(x), x, L)[0]
        
        for i in range(forcePointLoads.shape[1]):
            T += forcePointLoads[1,i] * forcePointLoads[2,i] * (1-np.heaviside(x-forcePointLoads[0,i], 1))
            
        for i in range(torquePointLoads.shape[1]):
            T += torquePointLoads[1,i] * (1-np.heaviside(x-torquePointLoads[0,i], 1))
        
        return T
        
    def plot(self, loading, lims, plots, subplots = True, step=0.01):
        # Ensure lims are of correct dimension
        assert len(lims) == 2
        
        x_vals = np.arange(lims[0], lims[1], step)
        
        axial_vals = []
        shear_vals = []
        moment_vals = []
        torsion_vals = []
        
        # ... 

# For testing
if __name__ == '__main__':
    calc = Calc(r'WP4\4.1\dataa0.txt', r'WP4\4.1\dataa10.txt')

    calc.lifttest()
    
    x_vals = np.arange(0, np.pi, 0.01)
    shear_vals = []
    moment_vals = []
    torsion_vals = []
    loading = lambda x: x
    lambda x: x**2
    
    
    for x in x_vals:
        shear_vals.append(calc.shear(x, np.pi, loading, NULL_ARRAY_2))
      
    for x in x_vals:  
        moment_vals.append(calc.moment(x, np.pi, calc.shear, NULL_ARRAY_2, 
                                       (np.pi, loading, NULL_ARRAY_2)))   
        
    for x in x_vals:  
        torsion_vals.append(calc.torsion(x, np.pi, lambda x: -1, lambda x: -1, lambda x: 0, NULL_ARRAY_3, NULL_ARRAY_2))
        
    plt.plot(x_vals, torsion_vals)
    plt.show()
