from calc import *
import numpy as np

if __name__ == '__main__':
    INT_INTRVLS = 100 # integration intervals
    calc = Calc(INT_INTRVLS)
    print(calc.calc_c())
    
    C_ROOT = 2*S/(b*1+b*TAPER_RATIO)
    print(C_ROOT, TAPER_RATIO*C_ROOT)
    print(b/2*np.tan(np.deg2rad(LAMBDA_LE)))
    