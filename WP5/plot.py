"""
Plotter for margin of safety.

This plotter takes an array of y values and corresponding applied stress values (alongside other parameters) to generate the plots.
Multiple plots can be generated in subplots by passing 2D arrays.

The parameters work as follows:
-yVals: numpy array, each row contains the y values for one subplot. Inhomogeneous arrays are allowed. Example:
        np.array([[y11, y12, y13], [y21, y22]])
-sigmaAppliedVals: numpy array, each row contains the aoplied stress values for one subplot. Inhomogeneous arrays are allowed. Example:
        np.array([[s11, s12, s13], [s21, s22]])
-sigma Fail: Failure stress, single number
-n: safety factor, single number
-dimSubplots: tuple, containing dimensions of how you want your subplots. This should match the input value arrays. Example:
        (2, 3) -> creates 2 rows and 3 columns of subplots
-titles: tuple, containing the title for each of your plots, in order.
-colors: tuple, containing the graph color for each of your plots, in order.
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

class DimensionError(Exception):
    pass


# y values, raw applied stress values, failure stress, safety factor for applied stress
def plotFailureMargin(yVals: NDArray, sigmaAppliedVals: NDArray, sigmaFail: float, n: float, dimSubplots: tuple, titles: tuple=(), colors: tuple=()) -> None:
    # Check for correct parameter dimensions/types
    if yVals.shape != sigmaAppliedVals.shape:
        raise DimensionError('yVals and sigmaAppliedVals must have same length!')
    if not isinstance(dimSubplots, tuple):
        raise TypeError('dimSubplots, titles must be tuples!')
        
    # Default values for colors
    if len(colors) == 0:
        colors = (int(yVals.shape[1]) * 'blue') 
        
    # To allow for inhomogeneous y/sigma arrays (plots with different # of data points),
    # numpy arrays are unpacked into python lists
    yArrays = [np.array(yVals[i]) for i in range(yVals.shape[0])]
    sigmaAppliedArrays = [np.array(yVals[i]) for i in range(yVals.shape[0])]
        
    # Array of failure stress values to allow for array division
    sigmaFailArrays = [n*np.full_like(sigmaAppliedArray, sigmaFail) for sigmaAppliedArray in sigmaAppliedArrays] 
    marginOfSafetyArrays = [sigmaFailArray/sigmaAppliedArray for sigmaFailArray, sigmaAppliedArray in zip(sigmaFailArrays, sigmaAppliedArrays)]
    
    # Plotting
    fig, axs = plt.subplots(*dimSubplots)
    
    for i, ax in enumerate(axs):
        ax.plot(yArrays[i], marginOfSafetyArrays[i], color=colors[i])
        ax.set_xlabel('y [m]')
        ax.set_ylabel('Margin of Safety [-]')
        ax.set_title(titles[i])
        ax.grid()
                
    fig.tight_layout()
    fig.suptitle('Margins of Safety')
    plt.show()
    

if __name__ == '__main__':
    plotFailureMargin(np.array([[1,2,3], [1,2]], dtype=object), 
                      np.array([[2,4,1], [3,5]], dtype=object), 
                      10,
                      1.5,
                      (1,2),
                      ('Doiuglaes', 'Test'),
                      ('red', 'blue'))
    
    
    
    
    
    
    
