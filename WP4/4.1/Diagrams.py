import matplotlib.pyplot as plt
import numpy as np

# Aircraft parameters

b_w = 20.4
L = b_w / 2
x = 0.0
dx = 0.1

# Moment diagram

x_vals = np.arange(0, L, 0.1)
y_vals = np.sin(x_vals)

# Graph
plt.plot(x_vals, y_vals)
plt.title('Moment Diagram')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)

plt.show()
