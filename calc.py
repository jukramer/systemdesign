import matplotlib.pyplot as plt
import numpy as np
from parameters import *

class Calc:
    def __init__(self):
        pass
    
    def WSMaxApproach(beta): # 1.23 gone for CS23
        return 1/beta*RHO_SL/2*(V_APP/1.23)**2*CL_MAX_L
        
    def WSMaxLField(beta):
        return 1/beta*LANDING_DIST/C_LFL*RHO_SL*CL_MAX_L/2
    
    def TSCruiseSpeed(beta,WS):
        pt = P_CR*(1+(GAMMA-1)/GAMMA*MACH_CRUISE**2)**(GAMMA/(GAMMA-1))   
        Tt = T_CR*(1+(GAMMA-1)/GAMMA*MACH_CRUISE**2)
        
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

        V_cr = MACH_CRUISE * np.sqrt(GAMMA*R*T_CR)

        return beta/alphat*((CD_0/2*RHO_CR*V_cr/(beta*WS))+(beta*WS/(np.pi*AR*e/2*RHO_CR*V_cr)))
        
           
    
            
    