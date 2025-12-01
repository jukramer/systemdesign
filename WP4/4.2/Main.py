import numpy as np
from Beam import *
np.set_printoptions(suppress=True)

points = [(0.2, 0.071507), (0.65, 0.071822), (0.65, -0.021653), (0.2, -0.034334)]

if True: # temporary
    y = np.linspace(0, 10, 100)
    def EllipticalLift(y, L, b):
        return (2*L/(np.pi*b)) * np.sqrt(np.maximum(0, 1 - (2*y/b)**2))

    def BendingMoment(y, L, b, n=2000):
        """
        Numerical bending moment at span station y using trapezoidal integration.
        y: spanwise station (0 to b/2)
        L: total lift
        b: total wingspan
        n: number of integration points
        """
        y_tip = b/2
        ys = np.linspace(y, y_tip, n)
        ls = EllipticalLift(ys, L, b)
        integrand = (ys - y) * ls
        return -np.trapezoid(integrand, ys)

    data = [(yi, BendingMoment(yi, L=7142*9.81, b=17.29)) for yi in y] # (y, M)
    plt.plot(y, [d[1]/np.abs(data[0][1]) for d in data], label=f'M | max {data[0][1]:.0f} Nm')
    data = np.array(data)

anything = Beam(100)
anything.load_wing_box(points=points, thickness=1/1000, root_chord=2.85, tip_chord=1.03, span=17.29)
anything.get_I_of_cross_section()
stringers = np.array([[0.05, 0.05**2], [0.015, 0.05**2]]) # [[z, A], [z, A]]
anything.get_displacement(data, E=73.1e9, stringers=stringers)
anything.get_twist(data, G=28e9)
anything.plot()

u a bitch micah