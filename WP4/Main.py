from WP4_1.calc import *
from WP4_2.Beam import *
import matplotlib.pyplot as plt
import numpy as np
import warnings

def main(debug=False):
    if not debug:
        warnings.simplefilter('ignore', category=UserWarning)
    
    # EXTERNAL LOADING
    calc = Calc(r'WP4\WP4_1\dataa0.txt', r'WP4\WP4_1\dataa10.txt')
    calc.set_load_case_from_flight(LOAD_FACTOR, W_MTOW)

    aeroLoading, inertialLoading, torsionLoading = lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[0], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[1], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[3]
    loadingDist = lambda x: calc.findLoadingDist(x)
    pointLoads, pointTorques = (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[2])(0), (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[4])(0)
    
    # INTERNAL LOADING
    y_data, M_data, T_data = calc.plot(aeroLoading,
                                       inertialLoading, 
                                       torsionLoading, 
                                       loadingDist, 
                                       pointLoads, 
                                       NULL_ARRAY_2, 
                                       pointTorques, 
                                       (0, HALF_SPAN),
                                       plot=True)
    
    return
    
    # DEFLECTION CALCULATIONS
    np.set_printoptions(suppress=True)

    plt.plot(y_data, M_data/np.abs(M_data[0]), label=f'M | max {M_data[0]:.0f} Nm')
    points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)]

    zis_is_ze_beam = Beam(y_data.size)
    zis_is_ze_beam.load_wing_box(points=points, thickness=1/1000, root_chord=2.85, tip_chord=1.03, span=17.29)
    zis_is_ze_beam.get_I_of_cross_section()
    stringers = np.array([[0.05, 0.05**2], [0.015, 0.05**2]])*0 # [[z, A], [z, A]]
    zis_is_ze_beam.get_displacement(np.column_stack((y_data, M_data)), E=72.4e9, stringers=stringers)
    zis_is_ze_beam.get_twist(np.column_stack((y_data, T_data)), G=28e9)
    zis_is_ze_beam.plot()


if __name__ == '__main__':
    main()