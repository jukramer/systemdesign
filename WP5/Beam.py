from globalParameters import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from Stringer import L_Stringer


class Beam():
    def __init__(self, stringers, intg_points: int = 100) -> None:
        self.intg_points = intg_points
        self.stringer_object: L_Stringer = stringers
        self.nBays = False
        
        # SPANWISE ARRAYS

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
        
        # self, NDArray, NDArray, 
        # SKIN WILL STAY CONST. THICKNESS!!
    def define_spanwise_arrays(self, y, posRibs, tStringersBay, bStringersBay, hStringersBay, nStringersBayTop, nStringersBayBottom):
        self.posRibs = np.array(posRibs)*HALF_SPAN # [% span] e.g. [0.1, 0.3, 0.5, ...] MUST Have 0 and 1 for tip and root!!!!
        assert int(posRibs[0]) == 0
        assert int(posRibs[-1]) == 1
        self.nBays = self.posRibs.shape[0] - 1
        sectionIDX = np.digitize(y, self.posRibs[1:-1])
        print(self.posRibs)
        
        self.distRibs = np.array([self.posRibs[i+1] - self.posRibs[i] for i in range(self.nBays)])[sectionIDX]
        self.nStringersTop = nStringersBayTop[sectionIDX] # TODO: THIS SHOULD WORK ?? TEST IT THO
        self.nStringersBottom = nStringersBayBottom[sectionIDX] # THIS SHOULD WORK ??
        self.tStringersBay = tStringersBay[sectionIDX]
        self.bStringersBay = bStringersBay[sectionIDX]
        self.hStringersBay = hStringersBay[sectionIDX]
        
        return self.distRibs, self.tStringersBay, self.bStringersBay, self.hStringersBay, self.nStringersTop, self.nStringersBottom
        
    def load_wing_box(self, points, stringer_count_top, stringer_count_bottom, aux_spar_endpoints, thickness, aux_spart_thickness, root_chord, tip_chord, span, posRibs=(0,)):
        self.points = points # [(x/c,z/c), ...] 
        self.aux_spar_endpoints = aux_spar_endpoints # [(x/c_start, y_start), (x/c_end, y_end)]
        self.thickness = thickness
        self.aux_spar_thickness = aux_spart_thickness
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.span = span

        self.define_stringers(self.points, stringer_count_top, stringer_count_bottom)

        # points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)] # [(x/c,z/c), ...] 
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
    
    def get_displacement(self, data, E, disable=False):
        s = self.intg_points
        y, M = data[:, 0], data[:, 1]
        c = self.get_chord(y)
        I = (self.Ixx_base_wingbox*c**3 # scaled wing box
            + self.stringer_object.Ixx
            + np.sum(self.stringer_object.area*(self.stringer_z_vals-self.centroid[1])**2, axis=0)*c**2 
            + np.where(y<=self.aux_spar_endpoints[1][1], 1/12*self.aux_spar_thickness*np.interp(y, [self.aux_spar_endpoints[0][1], self.aux_spar_endpoints[1][1]], [self.height_aux_spar_start, self.height_aux_spar_end])**3, 0)
            )
        self.Ixx_list = I
        
        d2v_dy2 = - M/(E * I)
        
        if disable:
            return
            
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

        # print(f'Volume: {self.volume:.4g} m³')
        return self.volume

    def report_stats(self):
        print(f'Deflected {self.v[-1]:.4g}m | Allowed {0.15*self.span:.4g}m')
        print(f'Twisted {self.theta[-1]*180/np.pi:.4g}° | Allowed {10.0:.4g}°')

    def konstantinos_konstantinopoulos(self, y, M, T=0, report=False):
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

    # TODO add torsion shear 
    def getShearStress(self, y, V):
        kV = 2 # shear factor, tau_max = kV*tau_avg (see reader app. F)
        # print(self.points)
        hFrontSpar = abs(self.points[0][1] - self.points[3][1])*self.get_chord(y)
        hRearSpar = abs(self.points[1][1] - self.points[2][1])*self.get_chord(y)
        
        return kV*V/(hFrontSpar*self.thickness+hRearSpar*self.thickness)
        
    # FAILURE STRESS CALCULATIONS
    def getFailureStresses(self, y):
        if not self.nBays:
            print('You must define the ribs first!')
            raise Exception
        
        # Failure stresses
        shearBuckStressCrit = self.shearBuckStress(y, )
        

        
        
         
    
    # Shear Buckling - this is a shear stress!!
    def shearBuckStress(self, y, rib_spacing):
        hFrontSpar = abs(self.points[0][1] - self.points[3][1])*self.get_chord(y)
        hRearSpar = abs(self.points[1][1] - self.points[2][1])*self.get_chord(y)
        heights = np.array([hFrontSpar, hRearSpar]).T

        a_b = rib_spacing/heights
        k_s = self.ShearBucklingInterpolation(a_b)

        t = self.thickness
        
        tau = np.pi**2 * k_s * E / (12*(1-POISSON_RATIO**2)) * (t/heights)**2
            
        return tau # shape(y.size, 2) for both spars
    
    def ShearBucklingInterpolation(self, a_b, plot=False):
        x_data=[1.00, 1.17, 1.50, 1.750, 2.00, 2.50, 3.0, 4.0, 5.0]
        y_data=[15.0, 13.0, 11.6, 10.84, 10.4, 9.84, 9.7, 9.5, 9.53]
        f=sp.interpolate.interp1d(x_data, y_data, kind='cubic', bounds_error=False, fill_value=(np.nan, 9.53))

        if plot:
            x_plt=np.arange(np.min(x_data)-1, np.max(x_data)+10, 0.001)
            y_plt=f(x_plt)
            img = plt.imread('WP5/old/Skin Buckling/web.png')  # replace with your image path
            fig, ax = plt.subplots()

            ax.imshow(img, extent=(-0.5, 5.55, 15.25, 3.85),
                origin='lower', aspect='auto', alpha=0.4, zorder=0)
            ax.plot(x_data, y_data, 'o', zorder=2, markersize=1, color='k')
            ax.plot(x_plt, y_plt, '-', zorder=3, linewidth=1)

            ax.set_xlim(-1, 6)
            ax.set_ylim(3, 17)
            plt.show()

        return f(a_b)
    
    # Skin buckling
    def skinBuckStress(self, y, rib_spacing):
        t = self.thickness
        chord = self.get_chord(y)
        edges = np.array(self.edge_lengths_list)
        b = np.outer(chord, edges)
        sigma = np.pi**2*self.SkinBucklingInterpolation(rib_spacing/b)*E / (12*(1-POISSON_RATIO**2)) * (t/b)**2
        return sigma
    
    def SkinBucklingInterpolation(self, a_b, plot=False):
        x_data=[0.7,  0.85, 1.0, 1.15, 1.3, 1.67, 2.0, 2.2, 2.5, 2.75, 3.0, 3.2, 3.35, 3.70, 4.18, 4.66, 4.8, 5.0]
        y_data=[10.7, 8.16, 6.8, 6.17, 5.8, 5.6,  4.9, 4.7, 4.6, 4.70, 4.5, 4.4, 4.4,  4.45, 4.27, 4.35, 4.3, 4.3]
        f=sp.interpolate.interp1d(x_data, y_data, kind='cubic', bounds_error=False, fill_value=(np.nan, 4.3))

        if plot:
            x_plt=np.arange(np.min(x_data)-1, np.max(x_data)+10, 0.001)
            y_plt=f(x_plt)
            img = plt.imread('WP5/old/Skin Buckling/image.png')  # replace with your image path
            fig, ax = plt.subplots()

            ax.imshow(img, extent=(-0.70, 6.0, 16.45, -1.05),
                origin='lower', aspect='auto', alpha=0.4, zorder=0)
            ax.plot(x_data, y_data, 'o', zorder=2, markersize=1, color='k')
            ax.plot(x_plt, y_plt, '-', zorder=3, linewidth=1)

            ax.set_xlim(-1, 16)
            ax.set_ylim(-1, 17)
            plt.show()

        return f(a_b)
    
    # Column Buckling - normal stress
    def colBuckStress(self, L):
        K = K_CC
        A = self.stringer_object.area
        I = self.Ixx_list
        return (K * np.pi**2 * E * I)/(L**2 * A)
    
    def calcStringerLen(self, sigma, K, I, A):
        return np.sqrt((K * np.pi**2 * E * I)/(sigma * A))
    
    def calcStringerArea(self, sigma, K, I, L):
        return np.sqrt((K * np.pi**2 * E * I)/(sigma * L**2))
    
    def calcStringerLenAll(self, sigma):
        # Wing tip / free end
        # TODO: I_stringer, A_Stringer
        if 1==1:
            raise RuntimeError('Check this')
        L_ribs_from_tip = self.calcStringerLen(sigma, K_FC, self.stringer_object.Ixx, self.stringer_object.area)
        # Between ribs / both fixed
        L_ribs_between = self.calcStringerLen(sigma, K_CC, self.stringer_object.Ixx, self.stringer_object.area)

        nRibs = np.ceil((b/2 - L_ribs_from_tip) / L_ribs_between).astype(int)

        return L_ribs_from_tip, L_ribs_between, nRibs
    
    def calcStringerAreaAll(self, sigma):
        LRibsFromTip, LRibsBetween, _ = self.calcStringerLenAll(sigma)

        if 1==1:
            raise RuntimeError('Check this')
        
        # Wing tip / free end
        A_ribs_from_tip = self.calcStringerArea(sigma, K_FC, self.stringer_object.Ixx, LRibsFromTip)

        # Between ribs / both fixed
        A_ribs_between = self.calcStringerArea(sigma, K_CC, self.stringer_object.Ixx, LRibsBetween)

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


if __name__=='__main__':
    wb = Beam(stringers=1, intg_points=865)
    # y = np.arange(0, 17.29/2, 17.29/2/wb.intg_points)
    # wb.shearBuckStress(y, 2)
    
    print(wb.define_spanwise_arrays(np.linspace(0, HALF_SPAN, 100),
                                    np.array([0, 0.5, 1]),
                                    np.array([1, 0.5]),
                                    np.array([1, 0.5]),
                                    np.array([1, 0.5]),
                                    np.array([2, 1]),
                                    np.array([2, 1]),))
    
