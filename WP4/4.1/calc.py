import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


class Calc():
    def __init__(self, file):
        with open(file) as data:
            dat = np.genfromtxt(data, skip_header=21, invalid_raise=False)

        ylst = dat[:,0]
        Cllst = dat[:,3]
        Cdlst = dat[:,5]
        Cmlst = dat[:,7]    

        print(ylst)
        print(Cllst)

        self.Cl = sp.interpolate.interp1d(ylst, Cllst, kind='cubic', fill_value='extrapolate')
        self.Cd = sp.interpolate.interp1d(ylst, Cdlst, kind='cubic', fill_value='extrapolate')
        self.Cm = sp.interpolate.interp1d(ylst, Cmlst, kind='cubic', fill_value='extrapolate')

        x_vals = np.arange(-7, 7, 0.05)

        # plt.plot(x_vals, self.Cl(x_vals))
        # plt.plot(x_vals, self.Cd(x_vals))
        # plt.plot(x_vals, self.Cm(x_vals))
        # plt.show()




    
    # Load Distribution as function of x
    def distrib(self):
        distrib = None

        return distrib
    
    # Normal force as function of x
    def axial(self, x, L, loading, pointLoads):
        N = sp.integrate.quad(loading, x, L)[0]
        
        for i in range(len(pointLoads)):
            N += pointLoads[i,1] * (1-np.heaviside(x-pointLoads[i,0], 1))
        
        return N

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
    

if __name__ == '__main__':
    calc = Calc(r'WP4\4.1\dataa0.txt')
    
    x_vals = np.arange(0, 10, 0.01)
    axial_vals = []
    loading = np.sin
    
    for x in x_vals:
        axial_vals.append(calc.axial(x, 10, loading, np.array([[2, 5],
                                                              [4, -2]])))
    
    # print(axial_vals := calc.axial(5, 10, loading, np.array([[2, 5],
    #                                                       [4, -2]])))
    
    plt.plot(x_vals, axial_vals)
    plt.show()    
    

    

