import numpy as np

class Beam():
    def __init__(self) -> None:
        pass

    def load_wing_box(self, points, thickness):
        self.points = points
        self.thickness = thickness

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
    
    def get_displacement(self, data, E):
        dvdz = [0,]

        for i in range(len(data)-1):
            y0, M0 = data(i)
            y1, M1 = data[min(i+1, len(data))]

            d_dvdz =  - (M1+M0)/2*(y1-y0) / (E * (self.Ixx*chord[i]**3 + self.Ixx*chord[i+1]**3)/2)


        

        