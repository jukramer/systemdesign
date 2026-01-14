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

order_of_mag = np.array([[1e-1, 1e-1, 1e-1, 1e0, 1e0, 1e2, 1e2]])

np.set_printoptions(suppress=True)

def x_from_print():
    # from optimisation run:
    x = np.array(
        [[0.001604, 0.004027, 0.002671, 0.03842, 0.03804, 19.86, 19.73],
        [0.0116, 0.00341, 0.002097, 0.04197, 0.04137, 19.35, 18.8],
        [0.05499, 0.003279, 0.001675, 0.04031, 0.04099, 19.38, 17.46],
        [0.107, 0.004129, 0.00118, 0.03246, 0.03049, 20, 18.34],
        [0.1476, 0.003667, 0.000745, 0.03095, 0.02923, 18.74, 16.65],
        [0.3293, 0.003324, 0.0005258, 0.0291, 0.02715, 18.34, 17.73],
        [0.4559, 0.003633, 0.0004261, 0.02776, 0.02615, 17.88, 17.48],
        [0.6238, 0.003086, 0.0003407, 0.02764, 0.02703, 17.91, 18.06],
        [0.8198, 0.003612, 0.000461, 0.02909, 0.02893, 18.71, 17.73]]
        )/order_of_mag

    # tweaked:
    x = np.array(
        [[0.000, 0.0042, 0.0011, 0.0105, 0.0105, 12., 12.],
        [0.0239, 0.0042, 0.0011, 0.0105, 0.0105, 12, 12.],
        [0.0616, 0.0042, 0.0011, 0.0105, 0.0105, 12., 12],
        [0.1122, 0.0042, 0.0011, 0.0105, 0.0105, 10., 10.],
        [0.1598, 0.0042, 0.0011, 0.0105, 0.0105, 10., 10.],
        [0.3235, 0.0042, 0.0011, 0.0105, 0.0105, 8.0, 8.0],
        [0.4562, 0.0042, 0.0011, 0.0105, 0.0105, 6.0, 6.0],
        [0.6231, 0.0042, 0.0011, 0.0105, 0.0105, 4., 4.],
        [0.8206, 0.0042, 0.0011, 0.0105, 0.0105, 3., 3.]]
        )/order_of_mag
    
    x = np.array(
        [[0.000, 0.0042, 0.0011, 0.0105, 0.0105, 12., 12.],
        [0.0500, 0.0042, 0.0011, 0.0105, 0.0105, 12., 12],
        [0.1122, 0.0042, 0.0011, 0.0105, 0.0105, 10., 10.],
        [0.1598, 0.0042, 0.0011, 0.0105, 0.0105, 10., 10.],
        [0.2700, 0.0042, 0.0011, 0.0105, 0.0105, 8.0, 8.0],
        [0.4250, 0.0042, 0.0011, 0.0105, 0.0105, 6.0, 6.0],
        [0.5900, 0.0042, 0.0011, 0.0105, 0.0105, 4.0, 4.0],
        [0.8206, 0.0042, 0.0011, 0.0105, 0.0105, 4.0, 4.0]]
        )/order_of_mag
    
    return x.T.flatten()

def prelim_x(idx):
# posRibs, tSkinBay, tStringersBay, bStringersBay, hStringersBay, nStringersBayTop, nStringersBayBottom = x
    x1 = np.array( 
        [[00, 1.6e-3, 1.0, 10.0, 10.0, 0, 0],
        [0.0, 1.6e-3, 1.0, 10.0, 10.0, 0, 0],
        [0.0, 1.6e-3, 1.0, 10.0, 10.0, 0, 0],
        [0.0, 1.6e-3, 1.0, 10.0, 10.0, 0, 0],
        [0.0, 1.6e-3, 1.0, 10.0, 10.0, 0, 0],
        [0.0, 1.6e-3, 1.0, 10.0, 10.0, 0, 0],
        [0.0, 1.6e-3, 1.0, 10.0, 10.0, 0, 0],
        [0.0, 1.6e-3, 1.0, 10.0, 10.0, 0, 0]]
        )/order_of_mag
    
    if idx==1:
        return x1
    
    x2 = np.array( 
        [[00, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2],
        [0.0, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2],
        [0.0, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2],
        [0.0, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2],
        [0.0, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2],
        [0.0, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2],
        [0.0, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2],
        [0.0, 1.e-3, 2.65e-3, 0.06, 0.06, 2, 2]]
        )/order_of_mag

    if idx==2:
        return x2
    
    x3 = np.array( 
        [[00, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13],
        [0.0, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13],
        [0.0, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13],
        [0.0, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13],
        [0.0, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13],
        [0.0, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13],
        [0.0, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13],
        [0.0, 1.e-3, 2.e-3, 0.012, 0.015, 13, 13]]
        )/order_of_mag

    if idx==3:
        return x3
    
    else:
        raise RuntimeError

def calc_mass(x):
    global bay_count, stressStack, marginArrayShear, marginArrayComp, marginArrayTens, modes, new_modes
    bay_count = x.size//7
    x = x.reshape(7, bay_count) * order_of_mag.T
    posRibs, tSkinBay, tStringersBay, bStringersBay, hStringersBay, nStringersBayTop, nStringersBayBottom = x
    
    posRibs = np.array(posRibs)
    tStringersBay = np.array(tStringersBay)
    tSkinBay = np.array(tSkinBay)
    bStringersBay = np.array(bStringersBay)
    hStringersBay = np.array(hStringersBay)
    nStringersBayTop = np.array(nStringersBayTop)
    nStringersBayBottom = np.array(nStringersBayBottom)
    
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 

    stringer_instance = L_Stringer(bStringersBay, hStringersBay, tStringersBay)
    # print(f'Stringer area: {stringer_instance.area[0]}')
    wb = Beam(stringers=stringer_instance, intg_points=865)
    wb.define_spanwise_arrays(y_data, posRibs, tStringersBay, bStringersBay, hStringersBay, nStringersBayTop, nStringersBayBottom)
    wb.load_wing_box(points=wing_box_points, thickness=tSkinBay, root_chord=2.85, tip_chord=1.03, span=17.29)
    v = wb.get_displacement(np.vstack((y_data, M_data)).T, E, False)
    theta = wb.get_twist(np.vstack((y_data, T_data)).T, G, False)
    wb.get_mass(y_data)
    
    # Applied Stresses
    global normalStressAppliedComp, normalStressAppliedTens, shearStressApplied

    normalStressAppliedTens = np.tile(np.maximum(0, np.max(wb.konstantinos_konstantinopoulos(y_data, M_data), axis=1, keepdims=True)), (1,2))
    normalStressAppliedComp = np.tile(np.minimum(0, np.min(wb.konstantinos_konstantinopoulos(y_data, M_data), axis=1, keepdims=True)), (1,6))
    shearStressApplied = wb.getShearStress(y_data, V_data, T_data)
    
    global stressStack
    # Critical Stresses
    stressStack = wb.getFailureStresses(y_data)[1]
    critStressArrayShear = stressStack[:, :2] # width 2
    critStressArrayComp = stressStack[:, 2:8] # width 6
    critStressArrayTens = stressStack[:, 8:] # width 2 

    modes = ['shearBuckStress', 'shearBuckStress', *['skinBuckStress' for _ in range(4)], 'colBuckStress', 'compYield', 'tensYield', 'crackPropStress']
    new_modes = ['Web shear buckling', 'Skin buckling', 'Stringer buckling', 'Compressive yield failure', 'Tensile yield failure', 'Crack propagation failure']
    
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

    rolled_ribs = np.roll(posRibs, -1)
    diff = rolled_ribs-posRibs
    diff[-1] = np.inf
    min_dist = np.min(diff)

    if iters % 100 == 0:
        print(f'Optimising ({iters})| {wb.mass:.1f}kg/{wb.volume:.3g}m³, Shear Margin: {minMarginShear:.3g}, Comp Margin: {minMarginComp:.3g}, Tens Margin: {minMarginTens:.3g}', end='\r')

    if iters % 500 == 0:
        arr = x.reshape(7, bay_count).T
        print('\n', np.array2string(arr, formatter={'float_kind':lambda v: f"{v:.4g}"}, separator=', '))
        wb.report_stats()


    # wb.plot()
    return wb.mass, np.array([minMarginShear, minMarginComp, minMarginTens, np.average([minMarginShear,minMarginComp,minMarginTens]), v[-1], theta[-1]*180/np.pi, min_dist])

def mass_wrap(x):
    global iters
    iters += 1
    # print(calcVol(x)[0])
    return calc_mass(x)[0]

def constr_wrap(x):
    # print(calcVol(x)[1].shape)
    return calc_mass(x)[1]

def optimise_main():
    global y_data, M_data, T_data, V_data, T_data
    global x_vals
    
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
    posRibs = (np.linspace(0, 1, 9)[:-1])**2
    ones = np.ones_like(posRibs)
    # initial_x = (np.array([posRibs, 4e-3*ones, 3e-3*ones, 5e-2*ones, 5e-2*ones, 20*ones, 20*ones])/order_of_mag.T).flatten()
    # x_vals = x_from_print()
    # print(initial_x.reshape(7, 9).T * order_of_mag)
    initial_x = (np.array([posRibs, 4e-3*ones, 2e-3*ones, 5e-2*ones, 5e-2*ones, np.linspace(40, 20, ones.size), np.linspace(40, 20, ones.size)])/order_of_mag.T).flatten()
    # print(initial_x.reshape(7, 8).T * order_of_mag)
    # quit()
    constraints_sigma = sp.optimize.NonlinearConstraint(constr_wrap, lb=[1, 1, 1, 1, -0.15*HALF_SPAN*2, -10.0, 0.01], ub=[np.inf, np.inf, np.inf, 1.5, 0.15*HALF_SPAN*2, 10.0, np.inf], keep_feasible=False)
    bounds_x = sp.optimize.Bounds((np.array([0*ones, 1e-3*ones,     1e-3*ones,     1e-2*ones,     1e-2*ones,     2*ones,   2*ones])/order_of_mag.T).flatten(), 
                                  (np.array([ones,   10e-3*ones, 10e-3*ones, 10e-2*ones, 10e-2*ones, 100*ones, 100*ones])/order_of_mag.T).flatten(),
                                  keep_feasible=False)
    
    optim = sp.optimize.minimize(mass_wrap, initial_x, method='trust-constr', bounds=bounds_x, constraints=constraints_sigma, ) # type:ignore | options = {"gtol": 1e-8}
    # optim = sp.optimize.shgo(volWrap, bounds=bounds_x, constraints=constraints_sigma) # type:ignore

    raw_result = optim.x
    print('')
    print("Success:", optim.success)
    print("Status:", optim.status)
    print("Message:", optim.message)
    print("Constraint violation:", optim.constr_violation)
    print("Optimality:", optim.optimality)
    print('\nResults:')
    arr = raw_result.reshape(7, bay_count).T * order_of_mag
    print(np.array2string(arr, formatter={'float_kind':lambda v: f"{v:.4g}"}, separator=', '))
    print(constr_wrap(raw_result))
    
    # print(stressList[0].shape)
    # print(stressList2[0].shape)
    # # plt.plot(iterList, stressList)
    # # plt.plot(iterList, stressList2)
    # # plt.yscale('log')
    # # plt.plot(np.linspace(0, 1, maxMarginArray.shape[0]), maxMarginArray)
    # plt.show()

def show_design():
    global y_data, M_data, T_data, V_data, T_data
    global x_vals
    
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
    x_vals = np.array( # Design one
        [[0.0000, 0.0042, 0.0011, 0.011, 0.011, 12, 12],
        [0.05000, 0.0042, 0.0011, 0.011, 0.011, 12, 12],
        [0.12000, 0.0042, 0.0011, 0.011, 0.011, 10, 10],
        [0.22672, 0.0042, 0.0011, 0.011, 0.011, 10, 10],
        [0.27000, 0.0042, 0.0011, 0.011, 0.011, 8, 8],
        [0.43000, 0.0042, 0.0011, 0.011, 0.011, 6, 6],
        [0.59500, 0.0042, 0.0011, 0.011, 0.011, 4, 4],
        [0.82000, 0.0042, 0.0011, 0.011, 0.011, 4, 4]]
        )/order_of_mag
    
    x_vals = np.array( # design two
        [[0.0000, 0.004, 0.0011, 0.024, 0.031, 35., 50.],
        [0.02800, 0.004, 0.0011, 0.024, 0.031, 35., 50.],
        [0.09000, 0.004, 0.0011, 0.024, 0.031, 35., 50.],
        [0.22672, 0.004, 0.0011, 0.024, 0.031, 35., 30.],
        [0.27000, 0.004, 0.0011, 0.024, 0.031, 20., 30.],
        [0.43350, 0.004, 0.0011, 0.024, 0.031, 20., 20.],
        [0.52830, 0.004, 0.0011, 0.024, 0.031, 20., 20.],
        [0.79260, 0.004, 0.0011, 0.024, 0.031, 14., 10.]]
        )/order_of_mag
    
    # x_vals = prelim_x(3)

    x_vals = x_vals.T.flatten()
    bay_count = x_vals.size//7
    # initial_x = (np.array([posRibs, 4e-3*ones, 2.7e-3*ones, 3e-2*ones, 3e-2*ones, np.linspace(20, 5, 9), np.linspace(20, 5, 9)])/order_of_mag.T).flatten()

    arr = x_vals.reshape(7, bay_count).T * order_of_mag
    # print(np.array2string(arr, formatter={'float_kind':lambda v: f"{v:.4g}"}, separator=', '))
    print(calc_mass(x_vals)[0])    
    plt.figure(figsize=(8,8))
    plt.plot(y_data/HALF_SPAN, marginArrayShear, color='red')
    plt.plot([19,20], [1,1], label=new_modes[0], color='red')
    plt.plot(y_data/HALF_SPAN, marginArrayComp[:, :4], color="#0000FF")
    plt.plot([19,20], [1,1], label=new_modes[1], color="#0000FF")

    if arr[0,-1] !=0:
        plt.plot(y_data/HALF_SPAN, marginArrayComp[:, 4], label=new_modes[2], color="#778EFF")
    plt.plot(y_data/HALF_SPAN, marginArrayComp[:, 5], label=new_modes[3], color="#77C9FF")
    plt.plot(y_data/HALF_SPAN, marginArrayTens[:, 0], color="#009E00")
    plt.plot([19,20], [1,1], label=new_modes[4], color="#009E00")
    plt.plot(y_data/HALF_SPAN, marginArrayTens[:, 1], color="#BCE052")
    plt.plot([19,20], [1,1], label=new_modes[5], color="#BCE052")

    plt.plot([0,20], [1,1], linestyle=':', color='k')
    plt.grid()
    plt.xlabel(r'$\frac{2y}{b}$ (fraction of half-span)')
    plt.ylabel(r'Allowed/applied stress (failure margin)')
    plt.xlim(0, 1)
    plt.ylim(0,5)
    plt.tight_layout()
    plt.legend()
    # plt.yscale('log')
    plt.show()

    # print(stressList[0].shape)
    # print(stressList2[0].shape)
    # # plt.plot(iterList, stressList)
    # # plt.plot(iterList, stressList2)
    # # plt.yscale('log')
    # # plt.plot(np.linspace(0, 1, maxMarginArray.shape[0]), maxMarginArray)
    # plt.show()

if __name__ == '__main__':
    # optimise_main() 
    show_design()