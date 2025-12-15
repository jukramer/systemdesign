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
    
    
    
    
    
    
    
