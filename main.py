from parameters import *
from calc import *

if __name__ == '__main__':
    # DRAWING PARAMS
    WS_MIN = 0
    WS_MAX = 4000 # program won't work if this is smaller than approach/landing constraints
    WS_INTRVL = 100
    
    calc = Calc()
    calc.drawMatchingDiagram(WS_MIN, WS_MAX, WS_INTRVL)
    # calc.designSpace(500, 100, 2601, 3599, None)