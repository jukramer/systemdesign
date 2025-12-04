from parameters import *
import matplotlib.pyplot as plt
import numpy as np

def plotDiagram(function):
    # Diagram

    z_vals = np.arange(0, HALF_SPAN, 0.1)
    y_vals = function(z_vals)
    line = np.zeros_like(z_vals)

    # Graph
    plt.plot(z_vals, y_vals, label='Moment', color='red')
    plt.plot(z_vals, line, linestyle='-', color='black')
    plt.title('Moment Diagram')
    plt.xlabel('z [m]')
    plt.ylabel('M [Nm]')
    plt.grid(True)
    plt.xlim(0,HALF_SPAN)

    plt.legend()
    plt.show()


