import matplotlib.pyplot as plt
from numpy.typing import NDArray

class DimensionError(Exception):
    pass

# y values, failure stress, raw applied stress values, safety factor for applied stress
def plotFailureMargin(yVals: NDArray, sigmaFail: float, sigmaApplied: NDArray, n: float):
    if yVals.shape != sigmaApplied.shape:
        raise DimensionError()
    
