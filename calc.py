import matplotlib.pyplot as plt
import numpy as np
from parameters import *

class Calc:
    def __init__(self):
        pass

    
    def alphaT(self, WS, T, P, rho, CL, MACH = None):
        V_alpha = np.sqrt(WS*2/rho*1/CL)
        
        if not MACH:
            MACH = V_alpha/np.sqrt(GAMMA*R*T)
        else:
            MACH = np.array(MACH)

        pt = P*(1+(GAMMA-1)/GAMMA*MACH**2)**(GAMMA/(GAMMA-1))
        Tt = T*(1+(GAMMA-1)/GAMMA*MACH**2)
        
        deltat = pt/P_SL
        thetat = Tt/T_SL
        alphat = None
    
        if B < 5:
            # if thetat < THETA_BREAK:
            #     alphat = deltat
            # elif thetat >= THETA_BREAK:
            #     alphat = deltat*(1-2.1*(thetat-THETA_BREAK)/THETA_BREAK)   
            alphat = np.where(thetat<THETA_BREAK, deltat, deltat*(1-2.1*(thetat-THETA_BREAK)/THETA_BREAK))
            
        elif B >= 5:
            # if theta < THETA_BREAK:
            #     alphat[i] = deltat[i]*(1-(0.43+0.014*B)*np.sqrt(MACH_CRUISE))
            # elif theta >= THETA_BREAK:
            #     alphat[i] = deltat[i]*(1-(0.43+0.014*B)*np.sqrt(MACH_CRUISE) - 3*(thetaval-THETA_BREAK)/(1.5+MACH_CRUISE))

            alphat = np.where(thetat < THETA_BREAK, deltat*(1-(0.43+0.014*B)*np.sqrt(MACH_CRUISE)), 
                              deltat*(1-(0.43+0.014*B)*np.sqrt(MACH_CRUISE) - 3*(thetat-THETA_BREAK)/(1.5+MACH_CRUISE)))


        return alphat
    
    def WSMaxApproach(self, beta): # /1.23 for FAR 25
        return 1/beta*RHO_SL/2*(V_APP/1.23)**2*CL_MAX_L
        
    def WSMaxLField(self, beta):
        return 1/beta*LANDING_DIST/C_LFL*RHO_SL*CL_MAX_L/2
    
    def TSCruiseSpeed(self, beta, WS):
        alphat = self.alphaT(WS, T_CR, P_CR, RHO_CR, CL_MAX_CR, MACH_CRUISE) # SL T/P Vals
        V_cr = MACH_CRUISE * np.sqrt(GAMMA*R*T_CR)
        print(V_cr)

        return beta/alphat*((CD_0/2*RHO_CR*V_cr**2/(beta*WS))+(beta*WS/(np.pi*AR*e/2*RHO_CR*V_cr**2)))

    def TSRateofClimb(self, beta, WS): # rate of climb
        alphat = self.alphaT(WS, T_SL, P_SL, RHO_SL, CL_MAX_TO) # SL T/P Vals
        return beta/alphat*(np.sqrt(c**2*RHO_SL/(2*beta*WS)*np.sqrt(CD_0*np.pi*AR*e))+2*np.sqrt(CD_0/(np.pi*AR*e)))

    def TSToF(self, beta, WS): # takeoff field length
        alphat = self.alphaT(WS, T_SL, P_SL, RHO_SL, CL_MAX_TO) # SL T/P Vals
        return 1.15*alphat*np.sqrt(WS/(TAKEOFF_DIST*KT*RHO_SL*9.80665*np.pi*AR*e)) + 4*H2/TAKEOFF_DIST

    def TSClimbGradient(self, beta, WS, CL, OEI, climbgrad):
        alphat = self.alphaT(WS, T_SL, P_SL, RHO_SL, CL) # SL T/P Vals
        
        V = np.sqrt(WS*2/RHO_SL/CL)
        
        if OEI:
            return N_E/(N_E-1)*beta/alphat*(climbgrad/V+2*np.sqrt(CD_0/(np.pi*AR*e)))
        else:
            return beta/alphat*(climbgrad/V+2*np.sqrt(CD_0/(np.pi*AR*e)))
    
    def drawMatchingDiagram(self, WS_MIN, WS_MAX):
        WS_VALS = np.arange(WS_MIN, WS_MAX, 100)
        
        # Approach/Landing Constraints
        WS_MAX_APP = self.WSMaxApproach(0.85) # 0.85 assumed
        WS_MAX_LAND = self.WSMaxLField(0.85)
        
        # TS Cruise/ROC/Takeoff
        TS_VALS_CR = self.TSCruiseSpeed(0.95, WS_VALS) # 0.95 assumed
        TS_VALS_ROC = self.TSRateofClimb(1, WS_VALS) # 1 assumed
        TS_VALS_TAKEOFF = self.TSToF(1, WS_VALS) # 1 assumed
        
        # TS Climb gradient
        TS_VALS_CLIMBGRAD_6 = self.TSClimbGradient(1, WS_VALS, CL_MAX_TO, True, 0)
        TS_VALS_CLIMBGRAD_7 = self.TSClimbGradient(1, WS_VALS, CL_MAX_TO, True, 0.024)
        TS_VALS_CLIMBGRAD_8 = self.TSClimbGradient(1, WS_VALS, CL_MAX_TO, True, 0.012)
        TS_VALS_CLIMBGRAD_9 = self.TSClimbGradient(1, WS_VALS, CL_MAX_L, False, 0.021)
        
        plt.xlim((WS_MIN, WS_MAX))
        
        plt.axvline(WS_MAX_APP, color='#236ee8', label='Stall', zorder=2)
        plt.axvline(WS_MAX_LAND, color="#be18f0", label='Landing Field', zorder=2)
        
        plt.plot(WS_VALS, TS_VALS_CR, color="#0FCBFF", label='Cruise', zorder=2)
        
        plt.plot(WS_VALS, TS_VALS_ROC, color="#1E6E10", label='Rate of Climb', zorder=2)
        plt.plot(WS_VALS, TS_VALS_TAKEOFF, color="#FF3C3C", label='Takeoff Field', zorder=2)
        plt.plot(WS_VALS, TS_VALS_CLIMBGRAD_6, color="#FF8A1D", label='Climb Gradient REQ_6', zorder=2)
        plt.plot(WS_VALS, TS_VALS_CLIMBGRAD_7, color="#FFE100", label='Climb Gradient REQ_7', zorder=2)
        plt.plot(WS_VALS, TS_VALS_CLIMBGRAD_8, color="#B1FF14", label='Climb Gradient REQ_8', zorder=2)
        plt.plot(WS_VALS, TS_VALS_CLIMBGRAD_9, color="#1FCF62", label='Climb Gradient REQ_9', zorder=2)
        
        a = np.maximum.reduce(np.array([TS_VALS_CR, TS_VALS_ROC, TS_VALS_TAKEOFF, TS_VALS_CLIMBGRAD_6, TS_VALS_CLIMBGRAD_7, TS_VALS_CLIMBGRAD_8, TS_VALS_CLIMBGRAD_9]))
        print(a)
        plt.fill_betweenx(a, 1, min(WS_MAX_APP, WS_MAX_LAND))
        
        plt.scatter(2084, 0.55, color='#ff0000', edgecolors='black', label='Design Point', zorder=3)
        
        plt.legend(title='Constraints')
        plt.xlabel(xlabel='Wing Loading [N/m$^2$]')
        plt.ylabel(ylabel='Thrust-to-Weight Ratio [-]')
    
        plt.show()
        
        # if I have time: function to draw hatchings next to lines