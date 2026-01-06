from Beam import *
from calc import *
from globalParameters import *
from Stringer import *
from plotSafetyMargin import plotFailureMargin
import matplotlib.pyplot as plt
import numpy as np
import warnings

np.set_printoptions(suppress=True)

def map_values(x):
    tString, tSkin, LStringBase, hString, nStringTop, nStringBottom, sRibs = x

    nStringTop_detach = int(nStringTop)
    nStringBottom_detach = int(nStringBottom)

    return tString, tSkin, LStringBase, hString, nStringTop_detach, nStringBottom_detach, sRibs

def calcVol(x):
    # tString, tSkin, LStringBase, hString, nStringTop, nStringBottom, sRibs = x
    # x[4] = initial_x[4]
    # x[5] = initial_x[5]
    tString_detach, tSkin_detach, LStringBase_detach, hString_detach, nStringTop_detach, nStringBottom_detach, sRibs_detach = map_values(x)
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 
    aux_spar_endpoints = [(0.425, -1), (0.2, -1)] # [(x/c_start, y_start), (x/c_end, y_end)] | Mind the units

    stringer_instance = L_Stringer(LStringBase_detach, hString_detach, tString_detach)
    wb = Beam(stringers=stringer_instance, intg_points=865)
    wb.load_wing_box(points=wing_box_points, stringer_count_top=nStringTop_detach, stringer_count_bottom=nStringBottom_detach, aux_spar_endpoints=aux_spar_endpoints, thickness=tSkin_detach, aux_spart_thickness=np.nan, root_chord=2.85, tip_chord=1.03, span=17.29)
    wb.get_displacement(np.vstack((y_data, M_data)).T, 1, True)
    vol = wb.get_volume()
    
    sigma_applied = wb.konstantinos_konstantinopoulos(y_data, M_data)
    sigma_skin = wb.skinBuckStress()
    sigma_col = np.repeat(wb.colBuckStress(sRibs_detach)[:, None], 4, 1)
    sigma_cr = np.minimum(sigma_skin, sigma_col)
    loss = np.sum(np.maximum(0, sigma_applied-sigma_cr))
    
    if loss > 0:
        print(f'Optimising stresses | {vol:.3f}m³, {loss:.1f}err          ', end='\r')
        # return loss*1e-5
    else:
        print(f'Optimising mass     | {vol:.3f}m³, {loss:.1f}err                 ', end='\r')
        # return vol
    
    return vol

def bucklingConstraints(x):
    tString_detach, tSkin_detach, LStringBase_detach, hString_detach, nStringTop_detach, nStringBottom_detach, sRibs_detach = map_values(x)
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 
    aux_spar_endpoints = [(0.425, -1), (0.2, -1)] # [(x/c_start, y_start), (x/c_end, y_end)] | Mind the units

    stringer_instance = L_Stringer(LStringBase_detach, hString_detach, tString_detach)
    wb = Beam(stringers=stringer_instance, intg_points=865)
    wb.load_wing_box(points=wing_box_points, stringer_count_top=nStringTop_detach, stringer_count_bottom=nStringBottom_detach, aux_spar_endpoints=aux_spar_endpoints, thickness=tSkin_detach, aux_spart_thickness=np.nan, root_chord=2.85, tip_chord=1.03, span=17.29)
    wb.get_displacement(np.vstack((y_data, M_data)).T, 1, True)
    vol = wb.get_volume()
    
    sigma_applied = wb.konstantinos_konstantinopoulos(y_data, M_data)/1e6
    sigma_skin = wb.skinBuckStress()/1e6
    sigma_col = np.repeat(wb.colBuckStress(sRibs_detach)[:, None], 4, 1)/1e6
    
    # print([np.sum(np.maximum(0, sigma_applied-sigma_skin)),
    #         np.sum(np.maximum(0, sigma_applied-sigma_col)),
    #         np.sum(np.maximum(0, sigma_applied-sigma_skin))], end='\r')
    
    return np.array([-np.sum(np.maximum(0, sigma_applied-sigma_skin)),
                     -np.sum(np.maximum(0, sigma_applied-sigma_col)),
                     -np.sum(np.maximum(0, sigma_applied-sigma_skin))])
    
# Wrapper functions for constraints (currently unused)
def bucklingConstraints1(x):
    return float(bucklingConstraints(x)[0])

def bucklingConstraints2(x):
    return float(bucklingConstraints(x)[1])

def bucklingConstraints3(x):
    return float(bucklingConstraints(x)[2])

def main(debug=False):   
    if not debug:
        warnings.simplefilter('ignore', category=UserWarning)

    # EXTERNAL LOADING
    calc = Calc(r'WP4\WP4_1\dataa0.txt', r'WP4\WP4_1\dataa10.txt')
    calc.set_load_case_from_flight(LOAD_FACTOR, W_MTOW)

    xVals = np.arange(0, HALF_SPAN, 0.01)
    
    aeroLoading, inertialLoading, torsionLoading = lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[0], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[1], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[3]
    loadingDist = lambda x: calc.findLoadingDist(x)
    pointLoads, pointTorques = (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[2])(0), (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[4])(0)
    
    # print(pointLoads)
    # print(pointTorques)
    
    # aeroLoadingVals = aeroLoading(xVals)
    # inertialLoadingVals = inertialLoading(xVals)
    # torsionLoadingVals = torsionLoading(xVals)
    
    # figLoad, (axLoad1, axLoad2, axLoad3) = plt.subplots(1,3)
    # axLoad1.plot(xVals, aeroLoadingVals, color='blue')
    # axLoad1.set_title('Distributed Normal Force')
    # axLoad1.set_xlabel('y [m]')
    # axLoad1.set_ylabel('Distributed Normal Force [N]')
    
    # axLoad2.plot(xVals, torsionLoadingVals, color='green')
    # axLoad2.set_title('Distributed Pitching Moment')
    # axLoad2.set_xlabel('y [m]')
    # axLoad2.set_ylabel('Distributed Pitching Moment [Nm]')
    
    # axLoad3.plot(xVals, inertialLoadingVals, color='red')
    # axLoad3.set_title('Distributed Inertial Loading')
    # axLoad3.set_xlabel('y [m]')
    # axLoad3.set_ylabel('Distributed Inertial Loading [N]')
    
    # figLoad.set_size_inches(15,5)
    # figLoad.suptitle(fr'External Loading', size='16', weight='semibold')
    # figLoad.tight_layout()
    # figLoad.savefig(fr'diagrams\Loading')
    # plt.show()

    # INTERNAL LOADING
    y_data, M_data, T_data = calc.plot(aeroLoading,
                                       inertialLoading, 
                                       torsionLoading, 
                                       loadingDist, 
                                       pointLoads, 
                                       NULL_ARRAY_2, 
                                       pointTorques, 
                                       (0, HALF_SPAN),
                                       subplots=True,
                                       plot=False)
        
    # DEFLECTION CALCULATIONS
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 

    aux_spar_endpoints = [(0.425, -1), (0.2, -1)] # [(x/c_start, y_start), (x/c_end, y_end)] | Mind the units
    # Thin stringers
    thin_stringer = L_Stringer(0.02, 0.02, 0.0015)
    stringer_count_top = 13
    stringer_count_bottom = 13
    skin_thickness = 1/1000

    # # Thick stringers
    # stringer_area = 0.312/1000 # m²
    # stringer_count_top = 2
    # stringer_count_bottom = 2
    # skin_thickness = 1/1000

    # # Thick skin
    # stringer_area = 0 # m²
    # stringer_count_top = 0
    # stringer_count_bottom = 0
    # skin_thickness = 1.6/1000

    many_stringer_beam = Beam(stringers=thin_stringer, intg_points=y_data.size)
    many_stringer_beam.load_wing_box(points=wing_box_points, stringer_count_top=stringer_count_top, stringer_count_bottom=stringer_count_bottom, aux_spar_endpoints=aux_spar_endpoints, thickness=skin_thickness, aux_spart_thickness=np.nan, root_chord=2.85, tip_chord=1.03, span=17.29)
    many_stringer_beam.get_displacement(np.column_stack((y_data, M_data)), E=72.4e9)
    many_stringer_beam.get_twist(np.column_stack((y_data, T_data)), G=28e9)
    many_stringer_beam.report_stats()
    many_stringer_beam.get_volume()
    many_stringer_beam.konstantinos_konstantinopoulos(y_data, T=T_data, M=M_data)
    many_stringer_beam.plot()

def optimise_main():
    global y_data, M_data
    global initial_x
    
    if not False:
        warnings.simplefilter('ignore', category=UserWarning)

    # EXTERNAL LOADING
    calc = Calc(r'WP4\WP4_1\dataa0.txt', r'WP4\WP4_1\dataa10.txt')
    calc.set_load_case_from_flight(LOAD_FACTOR, W_MTOW)

    xVals = np.arange(0, HALF_SPAN, 0.01)
    
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
                                       subplots=True,
                                       plot=False)
    
    initial_x = (5e-3, 5e-3, 0.1, 0.1, 13, 13, 2.0)
    bounds_x = [(0, 9e-3), (0, 9e-3), (0, 5e-2), (0, 5e-2), (0, 20), (0, 20), (0, 17)]
    constraints_sigma = [{'type': 'ineq', 'fun': bucklingConstraints1},
                         {'type': 'ineq', 'fun': bucklingConstraints2},
                         {'type': 'ineq', 'fun': bucklingConstraints3}]
    
    constraints_sigma = [{'type': 'ineq', 'fun': bucklingConstraints}]
    
    optim = sp.optimize.minimize(calcVol, initial_x, method='trust-constr', bounds=bounds_x, constraints=constraints_sigma)    
    raw_result = optim.x
    print('\nResult:')
    print(np.array(map_values(raw_result)))

if __name__ == '__main__':
    optimise_main()