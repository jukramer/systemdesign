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

        self.ylst0 = dat0[:,0]
        self.Cllst0 = dat0[:,3]
        Cdlst0 = dat0[:,5]
        Cmlst0 = dat0[:,7]
        clst0 = dat0[:,1]

        ylst10 = dat10[:,0]
        Cllst10 = dat10[:,3]
        Cdlst10 = dat10[:,5]
        Cmlst10 = dat10[:,7]
        clst10 = dat10[:,1]

        # Interpolate datapoints from XFLR5 data to obtain python functions

        self.Cl0 = sp.interpolate.interp1d(self.ylst0, self.Cllst0, kind='cubic', fill_value='extrapolate')
        self.Cd0 = sp.interpolate.interp1d(self.ylst0, Cdlst0, kind='cubic', fill_value='extrapolate')
        self.Cm0 = sp.interpolate.interp1d(self.ylst0, Cmlst0, kind='cubic', fill_value='extrapolate')

        self.Cl10 = sp.interpolate.interp1d(ylst10, Cllst10, kind='cubic', fill_value='extrapolate')
        self.Cd10 = sp.interpolate.interp1d(ylst10, Cdlst10, kind='cubic', fill_value='extrapolate')
        self.Cm10 = sp.interpolate.interp1d(ylst10, Cmlst10, kind='cubic', fill_value='extrapolate')

        # Global wing CL values and alpha
        self.CL0 = C_L0
        self.CL10 = C_L10
        self.alpha0 = 0.0
        self.alpha10 = 10.0

        # Placeholders for current load cases
        self.Cl = None
        self.Cd = None
        self.Cm = None
        self.alpha = None

        #self.lifttest()

    def set_load_case_CL(self, CLd):
        t = (CLd - self.CL0)/ (self.CL10 - self.CL0)

        self.Cl = lambda y: self.Cl0(y) + t * (self.Cl10(y) - self.Cl0(y))
        self.Cd = lambda y: self.Cd0(y) + t * (self.Cd10(y) - self.Cd0(y))
        self.Cm = lambda y: self.Cm0(y) + t * (self.Cm10(y) - self.Cm0(y))

        self.alpha = self.alpha0 + t * (self.alpha10 - self.alpha0)
        print(self.alpha)

    def set_load_case_from_flight(self, n, W, V=V_CR, rho=RHO, Sref=S):
        q_here = 0.5*rho*V**2
        CLd = n*W/(q_here*Sref)
        self.set_load_case_CL(CLd)

    def chord(self, y):
        c = C_ROOT + (C_TIP-C_ROOT) * 2 * y/b
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

    def calcNormal(self, alphaA):
        print('je moeder poept je moooooooooooooooooooooooeder - Larua')

    def inertialLoading(self, y, massWing, n=1):
        weightDens = n*g*massWing/S
        w = weightDens*self.chord(y)
        return w
    
    def pointLoading(self, thetaT, T):
        return T/2*np.sin(thetaT)
        
    ############ INTERNAL LOADING ##############

    # Normal force as function of x. pointLoads must have cols (position, load) (shape 2xn).
    # loading must be a python function.
    def axial(self, x, L, loading, pointLoads):
        N = sp.integrate.quad(loading, x, L)[0]

        for i in range(pointLoads.shape[1]):
            N += pointLoads[1,i] * (1-np.heaviside(x-pointLoads[0,i], 1))

        return N

    # Shear force as function of x. pointLoads must have cols (position, load) (shape 2xn)
    # loading must be a python function.
    def shear(self, x, L, loading, pointLoads):
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


    ############ PLOTTING ##############

    def plot(self, forceLoading, torsionLoading, loadingDist, pointLoads, pointMoments, pointTorques, lims, subplots = True, step=0.01):
        # Ensure arrays of correct dimension
        assert pointLoads.shape[0] == 3
        assert pointMoments.shape[0] == 2
        assert pointTorques.shape[0] == 2
        assert len(lims) == 2

        xMin, xMax = lims

        xVals = np.arange(xMin, xMax, step)

        shearVals = []
        momentVals = []
        torsionVals = []

        for x in xVals:
            shearVals.append(calc.shear(x, xMax, forceLoading, pointLoads))
            momentVals.append(calc.moment(x, xMax, calc.shear, pointMoments, (xMax, forceLoading, pointLoads)))
            torsionVals.append(calc.torsion(x, xMax, forceLoading, loadingDist, torsionLoading, pointLoads, pointTorques))

        # Plot with subplots
        if subplots:
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

            ax1.plot(xVals, shearVals)
            ax1.set_title('Shear Force Diagram')
            ax1.set_xlabel('y [m]')
            ax1.set_ylabel('Shear Force [m]')
            
            ax2.plot(xVals, momentVals)
            ax2.set_title('Bending Moment Diagram')
            ax2.set_xlabel('y [m]')
            ax2.set_ylabel('Bending Moment [Nm]')
            
            ax3.plot(xVals, torsionVals)
            ax3.set_title('Torsion Diagram')
            ax3.set_xlabel('y [m]')
            ax3.set_ylabel('Torsion [Nm]')
            
            plt.show()

        # Plot in sequential plots
        else:
            plt.plot(xVals, shearVals)
            plt.title('Shear Force Diagram')
            plt.xlabel('y [m]')
            plt.ylabel('Shear Force [m]')

            plt.show()
            plt.clf()

            plt.plot(xVals, momentVals)
            plt.title('Bending Moment Diagram')
            plt.xlabel('y [m]')
            plt.ylabel('Bending Moment [Nm]')

            plt.show()
            plt.clf()

            plt.plot(xVals, torsionVals)
            plt.title('Torsion Diagram')
            plt.xlabel('y [m]')
            plt.ylabel('Torsion [Nm]')

            plt.show()
            plt.clf()


if __name__ == '__main__':
    calc = Calc(r'WP4\4.1\dataa0.txt', r'WP4\4.1\dataa10.txt')

    xVals = np.arange(0, 8, 1e-2)
    clVals = calc.Cl0(xVals)

    # calc.lifttest()
    calc.set_load_case_from_flight(n_ult, W_MTOW)

    y = np.arange(0, HALF_SPAN, 1e-4)
    # L = calc.liftUnitSpan(y)
    # D = calc.dragUnitSpan(y)
    # M = calc.momentUnitSpan(y)

    # plt.plot(y, L)
    # plt.show()
    # plt.plot(y, D)
    # plt.show()
    # plt.plot(y, M)
    # plt.xlabel("y [m]")
    # plt.ylabel("Lift per unit span L'(y) [N/m]")
    # plt.grid(True)
    # plt.show()

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

    # plt.plot(x_vals, torsion_vals)
    # plt.show()
