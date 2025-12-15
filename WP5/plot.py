"""
Plotter for margin of safety.\n
\n
This plotter takes an array of y values and corresponding applied stress values (alongside other parameters) to generate the plots.\n
Multiple plots can be generated in subplots by passing 2D arrays.\n
\n
The parameters work as follows:\n
-yVals: numpy array, each row contains the y values for one subplot. Inhomogeneous arrays are allowed. Example:\n
        np.array([[y11, y12, y13], [y21, y22]])\n
-sigmaAppliedVals: numpy array, each row contains the aoplied stress values for one subplot. Inhomogeneous arrays are allowed. Example:\n
        np.array([[s11, s12, s13], [s21, s22]])\n
-sigma Fail: Failure stress, single number\n
-n: safety factor, single number\n
-dimSubplots: tuple, containing dimensions of how you want your subplots. This should match the input value arrays. Example:\n
        (2, 3) -> creates 2 rows and 3 columns of subplots\n
-titles: tuple, containing the title for each of your plots, in order (left-right then down).\n
-colors: tuple, containing the graph color for each of your plots, in order (left-right then down).
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
    if isinstance(axs[0,:], np.ndarray):
        axs = axs.ravel()
    
    for i, ax in enumerate(axs):
        # Handle unspecified colors/plots/titles
        try:
            ax.plot(yArrays[i], marginOfSafetyArrays[i], color=colors[i])
        except IndexError:
            try:
                ax.plot(yArrays[i], marginOfSafetyArrays[i], color='blue')
                print(f'Color for plot {i} not specified.')
                
            except IndexError:
                ax.set_axis_off()
        
        try:
            ax.set_title(titles[i])
        except IndexError:
            ax.set_title('')
            
        ax.set_xlabel('y [m]')
        ax.set_ylabel('Margin of Safety [-]')
        ax.grid()
                
    fig.tight_layout()
    fig.suptitle('Margins of Safety')
    plt.show()
    

if __name__ == '__main__':
    plotFailureMargin(np.array([[1,2,3], [1,2], [2,4,5,6,7]], dtype=object), 
                      np.array([[2,4,1], [3,5], [9,3,4,5,6]], dtype=object), 
                      10,
                      1.5,
                      (2,2),
                      ('Doiuglaes', 'Test'),
                      ('red', 'blue'))