import matplotlib.pyplot as plt
import numpy as np
from parameters import *
from calc import *

if __name__ == '__main__':
    calc = Calc()
    WSList = np.arange(100, 4000, 5)
    WSMaxApp = calc.WSMaxApproach(0.85)
    WSMaxFld = calc.WSMaxLField(0.85)
    TSCrsList = calc.TSCruiseSpeed(0.95, WSList)
