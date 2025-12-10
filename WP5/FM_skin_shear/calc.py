import numpy as np
import math
import os, sys
from parameters import *

current_folder = sys.path[0]
WP5 = os.path.dirname(os.path.dirname(current_folder))
sys.path.append(WP5)
from WP4.WP4_2 import Beam

def skin_shear_buckling(k_s, t, b):
   return np.pi**2 * k_s * E / (12*(1-poisson_ratio**2)) * (t/b)**2