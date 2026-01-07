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

order_of_mag = np.array([[1, 1e-3, 1e-3, 1e-2, 1e-2, 1e2, 1e2]])

np.set_printoptions(suppress=True)

def calcVol(x):
    bay_count = x.size//7
    x = x.reshape(7, bay_count) * order_of_mag.T
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
    v = wb.get_displacement(np.vstack((y_data, M_data)).T, E, True)
    theta = wb.get_twist(np.vstack((y_data, T_data)).T, G, True)

    vol = wb.get_volume()
    
    # Applied Stresses
    normalStressAppliedTens = np.tile(np.maximum(0, np.max(wb.konstantinos_konstantinopoulos(y_data, M_data), axis=1, keepdims=True)), (1,2))
    normalStressAppliedComp = np.tile(np.minimum(0, np.min(wb.konstantinos_konstantinopoulos(y_data, M_data), axis=1, keepdims=True)), (1,6))

    shearStressApplied = wb.getShearStress(y_data, V_data, T_data)
    
    # Critical Stresses
    stressStack = wb.getFailureStresses(y_data)[1]
    critStressArrayShear = stressStack[:, :2] # width 2
    critStressArrayComp = stressStack[:, 2:8] # width 6
    critStressArrayTens = stressStack[:, 8:] # width 2 

    modes = ['shearBuckStress', 'shearBuckStress', *['skinBuckStress' for _ in range(4)], 'colBuckStress', 'compYield', 'tensYield', 'crackPropStress']
    
    # Stress Margins
    marginArrayShear = critStressArrayShear/(shearStressApplied+1e-8)
    marginArrayComp = critStressArrayComp/(-normalStressAppliedComp+1e-8)
    marginArrayTens = critStressArrayTens/(normalStressAppliedTens+1e-8)

    minMarginShear = np.min(marginArrayShear)
    minMarginComp = np.min(marginArrayComp)
    argmin_comp = np.argmin(marginArrayComp, axis=1)[0]
    minMarginTens = np.min(marginArrayTens)
    
    # deltaArrayShear = critStressArrayShear - shearStressApplied
    # deltaArrayComp = critStressArrayComp - normalStressAppliedComp
    # deltaArrayTens = critStressArrayTens - normalStressAppliedTens
    
    global iters2
    iters2 += 1

    if iters2 % 50 == 0:
        print(f'Optimising ({iters2})| {vol:.3f}m³, Shear Margin: {minMarginShear}, Comp Margin: {minMarginComp}, {modes[argmin_comp]}, Tens Margin: {minMarginTens}', end='\r')

    if iters2 % 500 == 0:
        arr = x.reshape(7, bay_count).T
        print('\n', np.array2string(arr, formatter={'float_kind':lambda v: f"{v:.4g}"}, separator=', '))
    
    return vol, np.array([minMarginShear, minMarginComp, minMarginTens, v[-1], theta[-1]*180/np.pi])

def volWrap(x):
    # print(calcVol(x)[0])
    return calcVol(x)[0]

def constraintWrap(x):
    # print(calcVol(x)[1].shape)
    return calcVol(x)[1]

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
    posRibs = np.linspace(0, 1, 20)[:-1]
    ones = np.ones_like(posRibs)
    initial_x = (np.array([posRibs, 2e-3*ones, 3e-3*ones, 5e-2*ones, 5e-2*ones, 30*ones, 30*ones])/order_of_mag.T).flatten()
    
    constraints_sigma = sp.optimize.NonlinearConstraint(constraintWrap, lb=[1, 1, 1, -0.15*2*HALF_SPAN, -10.0], ub=[np.inf, np.inf, np.inf, 0.15*2*HALF_SPAN, 10.0])
    bounds_x = sp.optimize.Bounds((np.array([0*ones, 0*ones,     0*ones,     0*ones,     0*ones,     2*ones,   2*ones])/order_of_mag.T).flatten(), 
                                  (np.array([ones,   10e-3*ones, 10e-3*ones, 10e-2*ones, 10e-2*ones, 100*ones, 100*ones])/order_of_mag.T).flatten(),
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