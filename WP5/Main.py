from Beam import *
from calc import *
from globalParameters import *
from Stringer import *
from plotSafetyMargin import plotFailureMargin
import matplotlib.pyplot as plt
import numpy as np
import warnings

STRINGER_FACTOR = 1e5
iterList = []
stressList = []
stressList2 = []
iters = 1
iters2 = 1

np.set_printoptions(suppress=True)

def map_values(x):
    return x
    tString, tSkin, LStringBase, hString, nStringTop, nStringBottom, sRibs = x

    nStringTop_detach = int(nStringTop*STRINGER_FACTOR)
    nStringBottom_detach = int(nStringBottom*STRINGER_FACTOR)

    return tString, tSkin, LStringBase, hString, nStringTop_detach, nStringBottom_detach, sRibs

def calcVol(x):
    global maxMargin, maxMarginArray
    tString_detach, tSkin_detach, LStringBase_detach, hString_detach, sRibs_detach = map_values(x)
    nStringTop_detach, nStringBottom_detach = 20, 20
    # TODO: Stopped at calling function to define arrays
    
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 
    aux_spar_endpoints = [(0.425, -1), (0.2, -1)] # [(x/c_start, y_start), (x/c_end, y_end)] | Mind the units

    stringer_instance = L_Stringer(LStringBase_detach, hString_detach, tString_detach)
    wb = Beam(stringers=stringer_instance, intg_points=865)
    wb.load_wing_box(points=wing_box_points, thickness=tSkin_detach, root_chord=2.85, tip_chord=1.03, span=17.29)
    
    wb.get_displacement(np.vstack((y_data, M_data)).T, 1, True)
    vol = wb.get_volume()
    
    sigma_applied = wb.konstantinos_konstantinopoulos(y_data, M_data)
    sigma_skin = wb.skinBuckStress(y_data, sRibs_detach)
    # sigma_col = np.repeat(wb.colBuckStress(sRibs_detach)[:, None], 4, 1)
    sigma_col = sigma_skin
    sigma_cr = np.minimum(sigma_skin, sigma_col)
    
    maxMargin = np.max(np.maximum(sigma_applied/sigma_skin, sigma_applied/sigma_col))
    maxMarginArray = (np.maximum(sigma_applied/sigma_skin, sigma_applied/sigma_col))

    tau_web = wb.shearBuckStress(y_data, sRibs_detach)
    global iters2
    iters2+=1

    loss = np.sum(np.maximum(0, sigma_applied-sigma_cr))
    if iters2 % 100 == 0:
        if loss > 0:
            print(f'Optimising stresses | {vol:.3f}m³, {loss:.1f}err, {maxMargin:.2f}          ', end='\r')
        else:
            print(f'Optimising mass     | {vol:.3f}m³, {loss:.1f}err, {maxMargin:.2f}          ', end='\r')
    
    return vol*maxMargin 

def bucklingConstraints(x):
    global iters
    
    tString_detach, tSkin_detach, LStringBase_detach, hString_detach, sRibs_detach = map_values(x)
    nStringTop_detach, nStringBottom_detach = 20, 20
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 

    stringer_instance = L_Stringer(LStringBase_detach, hString_detach, tString_detach)
    wb = Beam(stringers=stringer_instance, intg_points=865)
    wb.load_wing_box(points=wing_box_points, thickness=tSkin_detach, root_chord=2.85, tip_chord=1.03, span=17.29)
    v = wb.get_displacement(np.vstack((y_data, M_data)).T, E, False)
    theta = wb.get_twist(np.vstack((y_data, T_data)).T, G)
    
    sigma_applied = wb.konstantinos_konstantinopoulos(y_data, M_data)
    sigma_skin = wb.skinBuckStress(y_data, sRibs_detach)
    sigma_col = np.repeat(wb.colBuckStress(sRibs_detach)[:, None], 4, 1)
    sigma_col=sigma_skin

    return np.array([np.max(sigma_applied - sigma_skin),
                     np.max(sigma_applied - sigma_col)])

def optimise_main():
    global y_data, M_data, T_data
    global initial_x
    global posRibs
    posRibs = (0,)
    
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
    
    initial_x = (4e-3, 2e-3, 9e-2, 9e-2, 2.0)
    constraints_sigma = sp.optimize.NonlinearConstraint(bucklingConstraints, lb=[-np.inf, -np.inf], ub=[0, 0])
    bounds_x = sp.optimize.Bounds([0, 0, 0, 0, 0], 
                                  [11e-3, 11e-3, 11e-2, 11e-2, 17],
                                  keep_feasible=True)
    
    optim = sp.optimize.minimize(calcVol, initial_x, method='trust-constr', bounds=bounds_x, constraints=constraints_sigma) # type:ignore   

    raw_result = optim.x
    print('')
    print("Success:", optim.success)
    print("Status:", optim.status)
    print("Message:", optim.message)
    print("Constraint violation:", optim.constr_violation)
    print("Optimality:", optim.optimality)
    print('\nResults:')
    print(np.array(map_values(raw_result)))
    
    print(stressList[0].shape)
    print(stressList2[0].shape)
    # plt.plot(iterList, stressList)
    # plt.plot(iterList, stressList2)
    # plt.yscale('log')
    plt.plot(np.linspace(0, 1, maxMarginArray.shape[0]), maxMarginArray)
    plt.show()

if __name__ == '__main__':
    optimise_main() 