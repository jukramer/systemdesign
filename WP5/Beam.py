from globalParameters import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


class Beam():
    def __init__(self, stringers, intg_points: int = 100) -> None:
        self.intg_points = intg_points
        self.stringer_object = stringers

    def define_stringers(self, wing_box_points, stringer_count_top, stringer_count_bottom):
        z_interp_top = lambda x: np.interp(x, [wing_box_points[0][0], wing_box_points[1][0]], [wing_box_points[0][1], wing_box_points[1][1]])
        z_interp_bottom = lambda x: np.interp(x, [wing_box_points[3][0], wing_box_points[2][0]], [wing_box_points[3][1], wing_box_points[2][1]])

        stringer_x_coords_top = [wing_box_points[0][0] + i/(stringer_count_top-1)*(wing_box_points[1][0]-wing_box_points[0][0]) for i in range(stringer_count_top)]
        stringer_x_coords_bottom = [wing_box_points[3][0] + i/(stringer_count_bottom-1)*(wing_box_points[2][0]-wing_box_points[3][0]) for i in range(stringer_count_bottom)]

        stringers_top = np.array([[z_interp_top(x), self.stringer_object.area] for x in stringer_x_coords_top]) if len(stringer_x_coords_top)>0 else np.array([[0,0]]) # [[z/c, A], [z/c, A]]
        stringers_bottom = np.array([[z_interp_bottom(x), self.stringer_object.area] for x in stringer_x_coords_bottom]) if len(stringer_x_coords_bottom)>0 else np.array([[0,0]]) # [[z/c, A], [z/c, A]]

        self.stringers = np.vstack((stringers_top, stringers_bottom)) # [[z/c, A], [z/c, A]]

        z_top = np.array([[z_interp_top(x)] for x in stringer_x_coords_top]) if len(stringer_x_coords_top)>0 else np.array([[0,0]])
        z_bottom = np.array([[z_interp_bottom(x)] for x in stringer_x_coords_bottom]) if len(stringer_x_coords_bottom)>0 else np.array([[0,0]])
        self.stringer_z_vals = np.vstack((z_top, z_bottom))

    def load_wing_box(self, points, stringer_count_top, stringer_count_bottom, aux_spar_endpoints, thickness, aux_spart_thickness, root_chord, tip_chord, span):
        self.points = points # [(x/c,z/c), ...] 
        self.aux_spar_endpoints = aux_spar_endpoints # [(x/c_start, y_start), (x/c_end, y_end)]
        self.thickness = thickness
        self.aux_spar_thickness = aux_spart_thickness
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.span = span

        self.define_stringers(self.points, stringer_count_top, stringer_count_bottom)

        # points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 
        # TODO: How are ribs/bays implemented ??
        chord_at_aux_spar_start = self.get_chord(aux_spar_endpoints[0][1])
        self.height_aux_spar_start = chord_at_aux_spar_start*(
            np.interp(aux_spar_endpoints[0][0], [points[0][0], points[1][0]], [points[0][1], points[1][1]]) # top
            - np.interp(aux_spar_endpoints[0][0], [points[3][0], points[2][0]], [points[3][1], points[2][1]])) # bottom
        
        chord_at_aux_spar_end = self.get_chord(aux_spar_endpoints[1][1])
        self.height_aux_spar_end = chord_at_aux_spar_end*(
            np.interp(aux_spar_endpoints[1][0], [points[0][0], points[1][0]], [points[0][1], points[1][1]]) # top
            - np.interp(aux_spar_endpoints[1][0], [points[3][0], points[2][0]], [points[3][1], points[2][1]])) # bottom
        
        self.get_I_of_cross_section()

    def get_I_of_cross_section(self):
        self.edge_centroids_list = []
        self.edge_lengths_list = []
        self.edge_angles_list = []

        for i in range(4): # for each edge
            x1, y1 = self.points[i % 4]
            x2, y2 = self.points[ (i+1) % 4]

            edge_centroid = np.array([(x2+x1)/2, (y2+y1)/2])
            edge_length = np.sqrt( (y2-y1)**2 + (x2-x1)**2 )
            edge_angle = np.arctan( (y2-y1)/(x2-x1) ) if abs(x2-x1)>1e-6 else np.pi/2

            self.edge_centroids_list.append(edge_centroid)
            self.edge_lengths_list.append(edge_length)
            self.edge_angles_list.append(edge_angle)

        self.centroid = np.array([0.0, 0.0])
        skin_area = sum(self.edge_lengths_list)*self.thickness
        stringer_area = np.sum(self.stringers[:,1])
        for i, c in enumerate(self.edge_centroids_list):
            self.centroid += c * self.edge_lengths_list[i]*self.thickness / (skin_area+stringer_area)
        for z, A in self.stringers:
            self.centroid[1] += z*A/(skin_area+stringer_area)

        self.Ixx_base_wingbox = 0
        self.Izz_base_wingbox = 0
        for i, c in enumerate(self.edge_centroids_list):
            self.Ixx_base_wingbox += (
                self.thickness * self.edge_lengths_list[i]**3 * np.sin(self.edge_angles_list[i])**2 / 12 # main component
                + self.edge_lengths_list[i] * self.thickness * (c-self.centroid)[1]**2                   # parallel axis component
                )

            self.Izz_base_wingbox += (
                self.thickness * self.edge_lengths_list[i]**3 * np.cos(self.edge_angles_list[i])**2 / 12 # main component
                + self.edge_lengths_list[i] * self.thickness * (c-self.centroid)[0]**2                   # parallel axis component
                )
            
        self.Ixz = 0

    def get_chord(self, y):
        frac = 2*np.abs(y)/self.span
        return (1-frac)*self.root_chord + (frac)*self.tip_chord
    
    def get_displacement(self, data, E):
        s = self.intg_points
        y, M = data[:, 0], data[:, 1]
        c = self.get_chord(y)
        I = (self.Ixx_base_wingbox*c**3 # scaled wing box
            + self.stringer_object.Ixx
            + np.sum(self.stringer_object.area*(self.stringer_z_vals-self.centroid[1])**2, axis=0)*c**2 
            + np.where(y<=self.aux_spar_endpoints[1][1], 1/12*self.aux_spar_thickness*np.interp(y, [self.aux_spar_endpoints[0][1], self.aux_spar_endpoints[1][1]], [self.height_aux_spar_start, self.height_aux_spar_end])**3, 0)
            )
        self.Ixx_list = I
        
        d2v_dy2 = - M / (E * I)

        integrand_1 = sp.interpolate.interp1d(y, d2v_dy2, kind='cubic', fill_value="extrapolate")
        dv_dy = np.empty(s)
        y2 = np.empty(s)
        for i in range(s):
            y2[i] = self.span/2*i/(s-1)
            dv_dy[i] = sp.integrate.quad(integrand_1, 0, self.span/2*i/(s-1))[0] # type: ignore

        integrand_2 = sp.interpolate.interp1d(y2, dv_dy, kind='cubic', fill_value="extrapolate")
        self.v = np.empty(s)
        for i in range(s):
            self.v[i] = sp.integrate.quad(integrand_2, 0, self.span/2*i/(dv_dy.size-1))[0] # type: ignore

        return self.v
    
    def get_twist(self, data, G):
        a, b, c, d = self.points # Order matters
        chord = self.get_chord(data[:, 0])
        self.Areas = (a[1]-d[1]+b[1]-c[1])/2*(c[0]-d[0]) * chord**2
        integral = chord*sum(self.edge_lengths_list)/self.thickness
        self.J = 4*self.Areas**2 / (integral)
        
        torque = data[:, 1]
        dtheta_dy = torque/(self.J*G)

        s = self.intg_points
        y = data[:, 0]
        integrand = sp.interpolate.interp1d(y, dtheta_dy, kind='cubic', fill_value="extrapolate")
        self.theta = np.empty(s)
        for i in range(s):
            self.theta[i] = sp.integrate.quad(integrand, 0, self.span/2*i/(s-1))[0] # type: ignore

        return self.theta

    def get_volume(self):
        y = np.linspace(0, self.span/2, self.intg_points)
        chord = self.get_chord(y)
        skin_area = sum(self.edge_lengths_list)*self.thickness*chord
        stringer_area = np.sum(self.stringers[:,1])

        integrand = sp.interpolate.interp1d(y, skin_area+stringer_area, kind='cubic', fill_value="extrapolate")
        self.volume = sp.integrate.quad(integrand, 0, self.span/2)[0] # type: ignore

        print(f'Volume: {self.volume:.4g} m³')
        return self.volume

    def report_stats(self):
        print(f'Deflected {self.v[-1]:.4g}m | Allowed {0.15*self.span:.4g}m')
        print(f'Twisted {self.theta[-1]*180/np.pi:.4g}° | Allowed {10.0:.4g}°')

    def konstantinos_konstantinopoulos(self, y, M, T, report=False):
        num_points = M.size
        self.normal_stress = np.empty((num_points, 4))
        chord = self.get_chord(y)

        for i, p in enumerate(self.points):
            x_c, z_c = p
            x = x_c * chord
            z = z_c * chord
            

            # self.normal_stress[:, i] = ((0*self.Ixx_list - M * self.Ixz)*x + (M*self.Izz - 0 * self.Ixz)*z) / (self.Ixx_list*self.Izz - self.Ixz**2)
            self.normal_stress[:, i] = M*z / self.Ixx_list
        if report:
            print(f'Max tensile: {np.max(self.normal_stress)/1e6:.0f}MPa, max compressive: {np.min(self.normal_stress)/1e6:.0f}MPa')

        return self.normal_stress

    # TODO: Function for shear stress

    # FAILURE STRESS CALCULATIONS
    # Shear Buckling - this is a shear stress!!
    def shearBuckStress(self, k_s, t, b):
        return np.pi**2 * k_s * E / (12*(1-POISSON_RATIO**2)) * (t/b)**2


    # Skin Buckling - normal stress
    def findkC(self,a,b): 
        x=[0.7,0.85,1.0,1.15,1.3,1.45,1.6,1.8,2.0,2.2,2.4,2.6,2.75,2.9,3.05,3.2,3.35,3.5,3.7,3.85,4.0,4.15,4.3,4.45,4.6,4.8,5.0]
        y=[10.7,7.5,6.8,6,5.8,5.5,5.5,5.1,4.9,4.7,4.5,4.5,4.6,4.5,4.4,4.4,4.4,4.3,4.3,4.3,4.2,4.2,4.2,4.1,4.1,4.1,4.1]
        f=sp.interpolate.interp1d(x,y,kind='cubic')
        xnew=np.arange(np.min(x),np.max(x),0.001)
        ynew=f(xnew)
        plt.plot(x,y,'o',xnew,ynew,'-')
        plt.show()
        return 0
    
    def skinBuckStress(self, t, b):
        return np.pi**2*self.findkC()*E / (12*(1-POISSON_RATIO**2)) * (t/b)**2
    
    # Column Buckling - normal stress
    def colBuckStress(self, K, A, L, I):
        return (K * np.pi**2 * E * I)/(L**2 * A)
    
    def calcStringerLen(self, sigma, K, I, A):
        return np.sqrt((K * np.pi**2 * E * I)/(sigma * A))
    
    def calcStringerArea(self, sigma, K, I, L):
        return np.sqrt((K * np.pi**2 * E * I)/(sigma * L**2))
    
    def calcStringerLenAll(self, sigma):
        # Wing tip / free end
        # TODO: I_stringer, A_Stringer
        L_ribs_from_tip = self.calcStringerLen(sigma, K_FC, I_Stringer, A_Stringer)
        # Between ribs / both fixed
        L_ribs_between = self.calcStringerLen(sigma, K_CC, I_Stringer, A_Stringer)

        nRibs = np.ceil((b/2 - L_ribs_from_tip) / L_ribs_between).astype(int)

        return L_ribs_from_tip, L_ribs_between, nRibs
    
    def calcStringerAreaAll(self, sigma):
        LRibsFromTip, LRibsBetween, _ = self.calcStringerLenAll()
        # Wing tip / free end
        A_ribs_from_tip = self.calcStringerArea(sigma, K_FC, I_Stringer, LRibsFromTip)

        # Between ribs / both fixed
        A_ribs_between = self.calcStringerArea(sigma, K_CC, I_Stringer, LRibsBetween)

        return A_ribs_from_tip, A_ribs_between

    # Crack Propagation
    def calcCCrit(self, sigma):
        return K_1C ** 2 / (np.pi * SHAPE_FACTOR ** 2 * sigma ** 2) * 10 ** 3 # [mm]

    # PLOTTING
    def plot(self):
        y = np.linspace(0, self.span/2, self.intg_points)

        plt.figure()
        plt.title('Stiffness Diagram')
        plt.plot(y, self.Ixx_list*1e4, label=f'I$_x$$_x$')
        plt.plot(y, self.J*1e4, label=f'J')
        plt.xlabel('y [$m$]')
        plt.ylabel('Stiffness × $10^4$ [$m^4$]')
        plt.xlim(0, self.span/2)
        plt.ylim(0, )
        plt.grid(which='both')
        plt.legend()
        plt.tight_layout()
        plt.show(block=False)
        
        plt.figure()
        plt.title('Deflection / Twisting Diagram')
        plt.plot(y, self.v/(np.max(np.abs(self.v))), label=f'v | max {self.v[-1]:.4g} m')
        plt.plot(y, self.theta/(np.max(np.abs(self.theta))), label=f'theta | max {np.max(np.abs(self.theta))*180/np.pi*np.sign(self.theta[-1]):.4g}°')
        plt.xlabel('y [$m$]')
        plt.ylabel('Normalised value')
        plt.xlim(0, self.span/2)
        plt.ylim(0, 1)
        plt.grid(which='both')
        plt.legend()
        plt.tight_layout()
        plt.show(block=True)