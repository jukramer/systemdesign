from parameters import *
from calc import *

if __name__ == '__main__':
    # DRAWING PARAMS
    WS_MIN = 100
    WS_MAX = 4000
    
    calc = Calc()
    calc.drawMatchingDiagram(WS_MIN, WS_MAX)