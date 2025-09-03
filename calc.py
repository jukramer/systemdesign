import matplotlib.pyplot as plt
import numpy as np
from parameters import *

class Calc:
    def __init__(self):
        pass
    
    def WSMaxApproach(self, beta): # 1.23 gone for CS23
        return 1/beta*RHO_SL/2*(V_APP/1.23)**2*CL_MAX_L
        
    def WSMaxLField(self, beta):
        return 1/beta*LANDING_DIST/C_LFL*RHO_SL*CL_MAX_L/2
    
    def TSCruiseSpeed(self, beta,  WS):
        pt = P_CR*(1+(GAMMA-1)/GAMMA*MACH_CRUISE**2)**(GAMMA/(GAMMA-1))   
        Tt = T_CR*(1+(GAMMA-1)/GAMMA*MACH_CRUISE**2)
        print(pt, Tt)
        deltat = pt/P_SL 
        thetat = Tt/T_SL
        alphat = None
        
        if B < 5:
            if thetat < THETA_BREAK:
                alphat = deltat
            elif thetat >= THETA_BREAK:
                alphat = deltat*(1-2.1*(thetat-THETA_BREAK)/THETA_BREAK)
        elif B >= 5:
            if thetat < THETA_BREAK:
                alphat = deltat*(1-(0.43+0.014*B)*np.sqrt(MACH_CRUISE))
            elif thetat >= THETA_BREAK:
                alphat = deltat*(1-(0.43+0.014*B)*np.sqrt(MACH_CRUISE) - 3*(thetat-THETA_BREAK)/(1.5+MACH_CRUISE))
        
        print(alphat)

        V_cr = MACH_CRUISE * np.sqrt(GAMMA*R*T_CR)
        print(V_cr)

        return beta/alphat*((CD_0/2*RHO_CR*V_cr**2/(beta*WS))+(beta*WS/(np.pi*AR*e/2*RHO_CR*V_cr**2)))


    def TSRateofClimb(self,beta,WS):
        return beta/alphat*(np.sqrt(c**2*RHO_SL/(2*beta*WS)*np.sqrt(CD_0*np.pi*AR*e))+2*np.sqrt(CD_0/(np.pi*AR*e)))



    
            
    