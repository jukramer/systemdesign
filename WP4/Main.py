import numpy as np
import os, sys
import matplotlib.pyplot as plt
from 4.2.Beam import *
np.set_printoptions(suppress=True)

WP1_data = np.load('WP4/4.2/Cases/case5.npz')
y_data = WP1_data['arr_0']
M_data = WP1_data['arr_1']
T_data = WP1_data['arr_2']
plt.plot(y_data, M_data/np.abs(M_data[0]), label=f'M | max {M_data[0]:.0f} Nm')
points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)]

zis_is_ze_beam = Beam(y_data.size)
zis_is_ze_beam.load_wing_box(points=points, thickness=1/1000, root_chord=2.85, tip_chord=1.03, span=17.29)
zis_is_ze_beam.get_I_of_cross_section()
stringers = np.array([[0.05, 0.05**2], [0.015, 0.05**2]])*0 # [[z, A], [z, A]]
zis_is_ze_beam.get_displacement(np.column_stack((y_data, M_data)), E=72.4e9, stringers=stringers)
zis_is_ze_beam.get_twist(np.column_stack((y_data, T_data)), G=28e9)
zis_is_ze_beam.plot()