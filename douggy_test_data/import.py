import numpy as np

modes = np.load('douggy_test_data/modes.npy')
StressStack = np.load('douggy_test_data/StressStack.npy')
y_vals = np.load('douggy_test_data/y_vals.npy')

print(modes, StressStack, y_vals)