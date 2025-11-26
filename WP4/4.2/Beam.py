import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate

if True: # dark mode
    plt.style.use('dark_background')
    plt.rcParams.update({
        'figure.facecolor': '#1a1a1a',   
        'axes.facecolor': '#1a1a1a',
        'axes.edgecolor': '#aaaaaa',      
        'axes.labelcolor': '#cccccc',
        'text.color': '#cccccc',
        'xtick.color': '#aaaaaa',
        'ytick.color': '#aaaaaa',
        'grid.color': '#444444',
        'lines.color': '#cccccc',
    })

class Beam():
    def __init__(self, intg_points: int = 100) -> None:
        self.v = []
        self.dvdz = []
        self.theta = []
        self.intg_points = intg_points

    def load_wing_box(self, points, thickness, root_chord, tip_chord, span):
        self.points = points
        self.thickness = thickness
        self.root_chord = root_chord
        self.tip_chord = tip_chord
        self.span = span

    def get_I_of_cross_section(self):
        self.edge_centroids_list = []
        self.edge_lengths_list = []
        self.edge_angles_list = []

        for i in range(4): # for each edge
            x1, y1 = self.points[i % 4]
            x2, y2 = self.points[ (i+1) % 4]

            edge_centroid = np.array([(x2+x1)/2, (y2+y1)/2])
            edge_length = np.sqrt( (y2-y1)**2 + (x2-x1)**2 )
            edge_angle = np.arctan( (y2-y1)/(x2-x1) ) if (x2-x1)>0 else np.pi/2

            self.edge_centroids_list.append(edge_centroid)
            self.edge_lengths_list.append(edge_length)
            self.edge_angles_list.append(edge_angle)

        self.centroid = 0
        for i, c in enumerate(self.edge_centroids_list):
            self.centroid += c * self.edge_lengths_list[i] / sum(self.edge_lengths_list)

        self.Ixx = 0
        self.Izz = 0
        for i, c in enumerate(self.edge_centroids_list):
            self.Ixx += (
                self.thickness * self.edge_lengths_list[i]**3 * np.sin(self.edge_angles_list[i])**2 / 12 # main component
                + self.edge_lengths_list[i] * self.thickness * (c-self.centroid)[1]**2                   # parallel axis component
                )

            self.Izz += (
                self.thickness * self.edge_lengths_list[i]**3 * np.cos(self.edge_angles_list[i])**2 / 12 # main component
                + self.edge_lengths_list[i] * self.thickness * (c-self.centroid)[0]**2                   # parallel axis component
                )
            
        self.Ixz = 0

        y = np.linspace(0.0, self.span/2, self.intg_points)
        chord = self.get_chord(y)
        self.Ixx_list = self.Ixx*chord**3

    def get_chord(self, y):
        frac = 2*np.abs(y)/self.span
        return (1-frac)*self.root_chord + (frac)*self.tip_chord
    
    def get_displacement(self, data, E, stringers):

        s = self.intg_points
        y, M = data[:, 0], data[:, 1]
        c = self.get_chord(y)
        d2v_dy2 = - M / (E * (self.Ixx*c**3 + np.sum(stringers[:, 1]*stringers[:, 0]**2, axis=0)*c**2 ))

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
        self.Areas = (a[1]-d[1]+b[1]-c[1])/2*(c[0]-d[0]) * self.get_chord(data[:, 0])**2
        chord = self.get_chord(data[:, 0])
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

    def plot(self):
        y = np.linspace(0, self.span/2, self.intg_points)
        plt.plot(y, self.v/(np.max(np.abs(self.v))), label=f'v | max {np.max(np.abs(self.v)):.4f} m')
        plt.plot(y, self.theta/(np.max(np.abs(self.theta))), label=f'theta | max {np.max(np.abs(self.theta))*180/np.pi*np.sign(self.theta[-1]):.4f} deg')
        plt.plot(y, self.Ixx_list/(np.max(np.abs(self.Ixx_list))), label=f'I$_x$$_x$ | max {np.max(np.abs(self.Ixx_list)):.4g} m$^4$')
        plt.plot(y, self.J/(np.max(np.abs(self.J))), label=f'J | max {np.max(np.abs(self.J)):.4g} m$^4$')
        plt.xlabel('y')
        plt.ylabel('value')
        plt.xlim(0, self.span/2)
        plt.ylim(-1, 1)
        plt.grid()
        plt.legend()
        plt.show()


        

        