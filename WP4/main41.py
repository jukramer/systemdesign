import numpy as np
import scipy as sp
from 4.1.calc import *

if __name__ == '__main__':
    calc = Calc(r'WP4\4.1\dataa0.txt', r'WP4\4.1\dataa10.txt')
                
    # External Loading    
    calc.set_load_case_from_flight(LOAD_FACTOR, W_MTOW)

    aeroLoading, inertialLoading, torsionLoading = lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[0], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[1], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[3]
    loadingDist = lambda x: calc.findLoadingDist(x)
    pointLoads, pointTorques = (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[2])(0), (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[4])(0)
    
    print(pointLoads)
    
    calc.plot(aeroLoading,
              inertialLoading, 
              torsionLoading, 
              loadingDist, 
              pointLoads, 
              NULL_ARRAY_2, 
              pointTorques, 
              (0, HALF_SPAN))
    
    
    
    
    

    