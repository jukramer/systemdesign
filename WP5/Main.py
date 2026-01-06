from Beam import *
from calc import *
from globalParameters import *
from Stringer import *
from plotSafetyMargin import plotFailureMargin
import matplotlib.pyplot as plt
import numpy as np
import warnings

class Optimizer:
    n1 = len(tStringVals := np.arange(1e-3, 4e-3, 4e-4))
    n2 = len(tSkinValues := np.arange(1e-3, 4e-3, 4e-4))
    n3 = len(LStringBaseVals := np.arange(5e-3, 10e-3, 5e-4))
    n4 = len(hStringVals := np.arange(5e-3, 10e-3, 5e-4))
    n5 = len(nStringTopVals := np.arange(4, 5, 1))
    n6 = len(nStringBotVals := np.arange(4, 5, 1))
    n7 = len(sRibVals := 33/np.arange(3, 6, 1))
    
    
    def __init__(self, y_data, M_data, T_data):
        self.y_data = y_data
        self.M_data = M_data
        self.T_data = T_data
        self.nTot = self.n1*self.n2*self.n3*self.n4*self.n5*self.n6*self.n7
        
    def setup_wing(self, tString, tSkin, LStringBase, hString, nStringTop, nStringBot, sRibs):
        wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 

        aux_spar_endpoints = [(0.425, -1), (0.2, -1)] # [(x/c_start, y_start), (x/c_end, y_end)] | Mind the units

        stringer_instance = L_Stringer(LStringBase, hString, tString)

        many_stringer_beam = Beam(stringers=stringer_instance, intg_points=self.y_data.size)
        many_stringer_beam.load_wing_box(points=wing_box_points, stringer_count_top=nStringTop, stringer_count_bottom=nStringBot, aux_spar_endpoints=aux_spar_endpoints, thickness=tSkin, aux_spart_thickness=np.nan, root_chord=2.85, tip_chord=1.03, span=17.29)
        # many_stringer_beam.get_displacement(np.column_stack((self.y_data, self.M_data)), E=72.4e9)
        # many_stringer_beam.get_twist(np.column_stack((self.y_data, self.T_data)), G=28e9)
        # many_stringer_beam.report_stats()
        many_stringer_beam.get_volume()
        # many_stringer_beam.konstantinos_konstantinopoulos(self.y_data, T=self.T_data, M=self.M_data)
        # many_stringer_beam.plot()
    
    def optimise(self):
        print(self.n1, self.n2, self.n3, self.n4, self.n5, self.n6, self.n7)
        i = 0
        
        for tString in self.tStringVals:
            for tSkin in self.tSkinValues:
                for LStringBase in self.LStringBaseVals:
                    for hString in self.hStringVals:
                        for nStringTop in self.nStringTopVals:
                            for nStringBottom in self.nStringBotVals:
                                for sRibs in self.sRibVals:
                                    i+=1
                                    print(f'{100*i/self.nTot}%')
                                    self.setup_wing(tString, tSkin, LStringBase, hString, nStringTop, nStringBottom, sRibs)

def calcVol(x):
    tString, tSkin, LStringBase, hString, nStringTop, nStringBottom, sRibs = x
    nStringTop = int(nStringTop)
    nStringBottom = int(nStringBottom)
    wing_box_points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 
    aux_spar_endpoints = [(0.425, -1), (0.2, -1)] # [(x/c_start, y_start), (x/c_end, y_end)] | Mind the units
    stringer_instance = L_Stringer(LStringBase, hString, tString)
    wb = Beam(stringers=stringer_instance, intg_points=865)
    wb.load_wing_box(points=wing_box_points, stringer_count_top=nStringTop, stringer_count_bottom=nStringBottom, aux_spar_endpoints=aux_spar_endpoints, thickness=tSkin, aux_spart_thickness=np.nan, root_chord=2.85, tip_chord=1.03, span=17.29)
    wb.get_displacement(np.vstack((y_data, M_data)).T, 1, True)
    mass = wb.get_volume()
    
    sigma = wb.konstantinos_konstantinopoulos(y_data, M_data)
    sigma_cr = wb.skinBuckStress()
    margin = np.sum(np.maximum(0, sigma-sigma_cr))*1e-9
    
    return mass + margin
    
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
    
    optim = Optimizer(y_data, M_data, T_data)
    optim = sp.optimize.minimize(calcVol, np.array([1e-3, 1e-3,1e-3,1e-3,5,5,5]))    
    print('Result:')
    print(optim.x)

if __name__ == '__main__':
    optimise_main()
    