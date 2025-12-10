from parameters import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

class Calc:
    def __init__(self):
        return

    def skinBucklingStress(t, b):
        return np.pi**2*kC*E / (12*(1-POISSION_RATIO**2)) * (t/b)**2


if __name__ == '__main__':
    pass