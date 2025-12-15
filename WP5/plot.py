import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

class DimensionError(Exception):
    pass


# y values, raw applied stress values, failure stress, safety factor for applied stress
def plotFailureMargin(yVals: NDArray, sigmaAppliedVals: NDArray, sigmaFail: float, n: float, dimSubplots: tuple, titles: tuple=(), colors: tuple=()) -> None:
    # Check for correct array dimensions
    if yVals.shape != sigmaAppliedVals.shape or not isinstance(dimSubplots, tuple):
        raise DimensionError('yVals and sigmaAppliedVals must have same length!')
        
    if len(colors) == 0:
        colors = (int(yVals.shape[1]) * 'blue') 
        
    sigmaFailArray = n*np.full_like(sigmaAppliedVals, sigmaFail)    
    marginOfSafety = sigmaFailArray/sigmaAppliedVals
    
    fig, axs = plt.subplots(*dimSubplots)
    
    for i, ax in enumerate(axs):
        print(yVals[i,:], marginOfSafety[i,:])
        ax.plot(yVals[i,:], marginOfSafety[i,:], color=colors[i])
        ax.set_xlabel('y [m]')
        ax.set_ylabel('Margin of Safety [-]')
        ax.set_title(titles[i])
        ax.grid()
                
    fig.tight_layout()
    fig.suptitle('Margins of Safety')
    plt.show()
    

if __name__ == '__main__':
    plotFailureMargin(np.array([[1,2,3], [1,2,3]]), 
                      np.array([[2,4,7], [3,5,6]]), 
                      10,
                      1.5,
                      (1,2),
                      ('Test', 'Test'),
                      ('red', 'blue'))
    
    
    
    
    
    
    
