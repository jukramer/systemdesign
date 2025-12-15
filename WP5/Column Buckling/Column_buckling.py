import numpy as np
from parameters import *

def sigma_crit(K, E, I, L, A):
    sigme_crit = (K * np.pi**2 * E * I)/(L**2 * A)
    return sigma_crit
