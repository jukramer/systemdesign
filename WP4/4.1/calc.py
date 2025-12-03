from parameters import *
from typing import Callable
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

        # Min and max x vals for plotting
        self.xMin, self.xMax = (0, 0)
        self.step = 1e-2
        
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

    def chord(self, y):
        return C_ROOT + (C_TIP-C_ROOT) * 2 * y/b

    def alpha_load_case(self, V, w, rho):
        L = w*n_ult
        C_L_d = L / (0.5*rho*(V**2)*S)

        return C_L_d

    def set_load_case_from_flight(self, n, W, V=V_CR, rho=RHO, Sref=S):
        q_here = 0.5*rho*V**2
        CLd = n*W/(q_here*Sref)
        
        t = (CLd - self.CL0)/ (self.CL10 - self.CL0)

        self.Cl = lambda y: self.Cl0(y) + t * (self.Cl10(y) - self.Cl0(y))
        self.Cd = lambda y: self.Cd0(y) + t * (self.Cd10(y) - self.Cd0(y))
        self.Cm = lambda y: self.Cm0(y) + t * (self.Cm10(y) - self.Cm0(y))

        self.alpha = self.alpha0 + t * (self.alpha10 - self.alpha0)

    ######### EXTERNAL LOADING ##############
    # AERODYNAMIC LOADING
    def calcNormal(self, y, alphaA):
        c = self.chord(y)
        alphaA = np.deg2rad(alphaA)
        L = self.Cl(y)*q*c
        D = self.Cd(y)*q*c
        
        return np.cos(alphaA)*L + np.sin(alphaA)*D
    
    def momentUnitSpan(self, y):
        return self.Cm(y)*q*self.chord(y)**2

    # INERTIAL LOADING
    def inertialLoading(self, y, massWing, n=1):
        weightDens = n*g*massWing/S
        return -weightDens*self.chord(y)
        
    # PROPULSIVE LOADING
    def propulsiveLoading(self, thetaT, T):
        return np.array(([0], [T/2*np.sin(thetaT)], [0.25*C_ROOT]))
    
    def propulsiveMoment(self, thetaT, T, d):
        return np.array(([0], [T/2*np.sin(thetaT)*d]))
    
    # TOTAL LOADING
    def totalLoading(self, x, n, mWing):
        return self.calcNormal(x, self.alpha-WING_TRIM) + self.inertialLoading(x, mWing, n), self.propulsiveLoading(self.alpha - WING_TRIM, T_TO), self.momentUnitSpan(x), self.propulsiveMoment(self.alpha - WING_TRIM, T_TO, d_prop)
        
    ############ INTERNAL LOADING ##############
    # Shear force as function of x. pointLoads must have cols (position, load) (shape 2xn)
    # loading must be a python function.
    def shear(self, x, loading: Callable, pointLoads):
        xValsInt = np.arange(x, self.xMax, self.step)
        yValsInt = loading(xValsInt)
        
        V = sp.integrate.simpson(x=xValsInt, y=yValsInt)

        for i in range(pointLoads.shape[1]):
            V += pointLoads[1,i] * (1-np.heaviside(x-pointLoads[0,i], 1))

        return V

    # Moment force as function of x. pointLoads must have cols (position, load) (shape 2xn)
    # loading must be a python function.
    def moment(self, x, shearVals, pointMoments):
        xValsInt = np.arange(x, self.xMax, self.step)
        i = shearVals.shape[0] - xValsInt.shape[0]
        M = sp.integrate.simpson(x=xValsInt, y=shearVals[i:])

        for i in range(pointMoments.shape[1]):
            M += pointMoments[1,i] * (1-np.heaviside(x-pointMoments[0,i], 1))

        return -M

    # Torsion as function of x. forcePointLoads must have cols (position, load, dist) (shape 3xn);
    # torquePointLoads must have cols (position, load) (shape 2xn)
    # forceloading, torsionLoading, and loadingDist must be a python function.
    def torsion(self, x, loadingVals, forcePointLoads, torquePointLoads):
        xValsInt = np.arange(x, self.xMax, self.step)
        i = loadingVals.shape[0] - xValsInt.shape[0]
        T = sp.integrate.simpson(x=xValsInt, y=loadingVals[i:])
        
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

        self.xMin, self.xMax = lims
        self.step = step
        xVals = np.arange(self.xMin, self.xMax, step)        

        self.shearVec = np.vectorize(self.shear, signature='(),(),(3,1)->()')
        shearVals = self.shearVec(xVals, forceLoading, pointLoads)
        
        self.momentVec = np.vectorize(self.moment, signature='(),(n),(m,l)->()')
        momentVals = self.momentVec(xVals, shearVals, pointMoments)
        
        self.torsionVec = np.vectorize(self.torsion, signature='(),(m),(3,1),(2,1)->()')
        torsionLoadVals = torsionLoading(xVals) + forceLoading(xVals)*loadingDist(xVals)
        torsionVals = self.torsionVec(xVals, torsionLoadVals, pointLoads, pointTorques)
        
        np.savez('shearVals', x=xVals, y=shearVals)
        np.savez('momentvals', x=xVals, y=momentVals)
        np.savez('torsionVals', x=xVals, y=torsionVals)

        print('Plotting!')
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

            # plt.plot(xVals, momentVals)
            # plt.title('Bending Moment Diagram')
            # plt.xlabel('y [m]')
            # plt.ylabel('Bending Moment [Nm]')

            # plt.show()
            # plt.clf()

            # plt.plot(xVals, torsionVals)
            # plt.title('Torsion Diagram')
            # plt.xlabel('y [m]')
            # plt.ylabel('Torsion [Nm]')

            # plt.show()
            # plt.clf()


if __name__ == '__main__':
    calc = Calc(r'WP4\4.1\dataa0.txt', r'WP4\4.1\dataa10.txt')

    calc.set_load_case_from_flight(n_ult, W_MTOW)
        

    wlst = [W_MTOW, W_minusfuel, W_OEM]
    Vlst = [1.5*V_CR, V_CR, V_stallwflaps]
    RHOlst = [RHO, RHO_SL]
    for V in Vlst:
        for w in wlst:
            for rho in RHOlst:
                CLd = calc.alpha_load_case(V, w, rho)
                print(
                f"V = {V:6.2f} m/s | "
                f"W = {w:8.0f} N ({w/g:6.1f} kg) | "
                f"rho = {rho:5.3f} kg/m³ | "
                f"CLd = {CLd:8.3f}"
                )

    forceLoading, torsionLoading = lambda x: calc.totalLoading(x, 1, M_WING)[0], lambda x: calc.totalLoading(x, 1, 1000)[2]
    loadingDist = lambda x: calc.chord(x)
    pointLoads, pointTorques = (lambda x: calc.totalLoading(x, 1, M_WING)[1])(0), (lambda x: calc.totalLoading(x, 1, 1000)[3])(0)
    
    print(pointLoads)
    
    calc.plot(forceLoading, 
              torsionLoading, 
              loadingDist, 
              pointLoads, 
              NULL_ARRAY_2, 
              pointTorques, 
              (0, HALF_SPAN))
    
    