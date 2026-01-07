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

def calcVol(x):
    posRibs, tStringersBay, tSkinBay, bStringersBay, hStringersBay, nStringersBayTop, nStringersBayBottom = x
    
    posRibs = np.array(posRibs)
    tStringersBay = np.array(tStringersBay)
    tSkinBay = np.array(tSkinBay)
    bStringersBay = np.array(bStringersBay)
    hStringersBay = np.array(hStringersBay)
    nStringersBayTop = np.array(nStringersBayTop)
    nStringersBayBottom = np.array(nStringersBayBottom)
    
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 

    stringer_instance = L_Stringer(bStringersBay, hStringersBay, tStringersBay)
    wb = Beam(stringers=stringer_instance, intg_points=865)
    wb.define_spanwise_arrays(y_data, posRibs, tStringersBay, bStringersBay, hStringersBay, nStringersBayTop, nStringersBayBottom)
    wb.load_wing_box(points=wing_box_points, thickness=tSkinBay, root_chord=2.85, tip_chord=1.03, span=17.29)
    v = wb.get_displacement(np.vstack((y_data, M_data)).T, E, False)
    theta = wb.get_twist(np.vstack((y_data, T_data)).T, G)

    vol = wb.get_volume()
    
    # Applied Stresses
    normalStressAppliedTens = np.tile(np.maximum(0, np.max(wb.konstantinos_konstantinopoulos(y_data, M_data), axis=1)), (1,2))
    normalStressAppliedComp = np.tile(np.minimum(0, np.min(wb.konstantinos_konstantinopoulos(y_data, M_data), axis=1)), (1,3))
    shearStressApplied = wb.getShearStress(y_data, V_data, T_data)
    
    # Critical Stresses
    critStressArrayShear = wb.getFailureStresses(y_data)[1][0,:] # width 1
    critStressArrayComp = wb.getFailureStresses(y_data)[1][1:4,:] # width 3
    critStressArrayTens = wb.getFailureStresses(y_data)[1][4:,:] # width 2 
    
    # Stress Margins
    marginArrayShear = critStressArrayShear/shearStressApplied
    marginArrayComp = critStressArrayComp/normalStressAppliedComp
    marginArrayTens = critStressArrayTens/normalStressAppliedTens
    
    maxMarginShear = np.max(marginArrayShear)
    maxMarginComp = np.max(marginArrayComp)
    maxMarginTens = np.max(marginArrayTens)
    
    # deltaArrayShear = critStressArrayShear - shearStressApplied
    # deltaArrayComp = critStressArrayComp - normalStressAppliedComp
    # deltaArrayTens = critStressArrayTens - normalStressAppliedTens
    
    global iters2
    iters2 += 1

    if iters2 % 100 == 0:
        print(f'Optimising | {vol:.3f}m³, Shear Margin: {maxMarginShear}, Comp Margin: {maxMarginComp}, Tens Margin: {maxMarginTens}', end='\r')
    
    return vol, np.array([maxMarginShear, maxMarginComp, maxMarginTens, v[-1], theta[-1]*180/np.pi])

def volWrap(x):
    print(calcVol(x)[0].shape)
    return calcVol(x)[0]

def constraintWrap(x):
    print(calcVol(x)[1].shape)
    return calcVol(x)[1]

# def bucklingConstraints(x):
#     global iters
    
#     tString_detach, tSkin_detach, LStringBase_detach, hString_detach, sRibs_detach = map_values(x)
#     nStringTop_detach, nStringBottom_detach = 20, 20
#     wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 

#     stringer_instance = L_Stringer(LStringBase_detach, hString_detach, tString_detach)
#     wb = Beam(stringers=stringer_instance, intg_points=865)
#     wb.define_spanwise_arrays()
#     wb.load_wing_box(points=wing_box_points, thickness=tSkin_detach, root_chord=2.85, tip_chord=1.03, span=17.29)
#     v = wb.get_displacement(np.vstack((y_data, M_data)).T, E, False)
#     theta = wb.get_twist(np.vstack((y_data, T_data)).T, G)
    
#     sigma_applied = wb.konstantinos_konstantinopoulos(y_data, M_data)
#     sigma_skin = wb.skinBuckStress(y_data, sRibs_detach)
#     sigma_col = np.repeat(wb.colBuckStress(sRibs_detach)[:, None], 4, 1)
#     sigma_col=sigma_skin

#     return np.array([np.max(sigma_applied - sigma_skin),
#                      np.max(sigma_applied - sigma_col),
#                      v[-1], # type:ignore
#                      theta[-1]*180/np.pi])

def optimise_main():
    global y_data, M_data, T_data, V_data, T_data
    global initial_x
    
    if not False:
        warnings.simplefilter('ignore', category=UserWarning)

    # EXTERNAL LOADING
    calc = Calc(r'WP4\WP4_1\dataa0.txt', r'WP4\WP4_1\dataa10.txt')
    calc.set_load_case_from_flight(LOAD_FACTOR, W_MTOW)

    aeroLoading, inertialLoading, torsionLoading = lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[0], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[1], lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[3]
    loadingDist = lambda x: calc.findLoadingDist(x)
    pointLoads, pointTorques = (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[2])(0), (lambda x: calc.totalLoading(x, LOAD_FACTOR, M_WING)[4])(0)
    
    # INTERNAL LOADING
    y_data, M_data, T_data, V_data = calc.plot(aeroLoading,
                                       inertialLoading, 
                                       torsionLoading, 
                                       loadingDist, 
                                       pointLoads, 
                                       NULL_ARRAY_2, 
                                       pointTorques, 
                                       (0, HALF_SPAN),
                                       subplots=True,
                                       plot=False)
    
    # OPTIMIZATION
    # posRibs, tSkinBay, tStringersBay, bStringersBay, hStringersBay, nStringersBayTop, nStringersBayBottom = x
    initial_x = (np.array([0,1]), 1e-3, 1e-3, 5e-2, 5e-2, 10, 10) 
    
    constraints_sigma = sp.optimize.NonlinearConstraint(constraintWrap, lb=[1, 1, 1, -0.15*2*HALF_SPAN, -10.0], ub=[np.inf, np.inf, np.inf, 0.15*2*HALF_SPAN, 10.0])
    bounds_x = sp.optimize.Bounds([0,         0,     0,     0,     0,     2,  2], 
                                  [HALF_SPAN, 10e-3, 10e-3, 10e-2, 10e-2, 17, 17],
                                  keep_feasible=True)
    
    optim = sp.optimize.minimize(volWrap, initial_x, method='trust-constr', bounds=bounds_x, constraints=constraints_sigma) # type:ignore   

    raw_result = optim.x
    print('')
    print("Success:", optim.success)
    print("Status:", optim.status)
    print("Message:", optim.message)
    print("Constraint violation:", optim.constr_violation)
    print("Optimality:", optim.optimality)
    print('\nResults:')
    print(np.array(raw_result))
    
    # print(stressList[0].shape)
    # print(stressList2[0].shape)
    # # plt.plot(iterList, stressList)
    # # plt.plot(iterList, stressList2)
    # # plt.yscale('log')
    # # plt.plot(np.linspace(0, 1, maxMarginArray.shape[0]), maxMarginArray)
    # plt.show()

if __name__ == '__main__':
    optimise_main() 