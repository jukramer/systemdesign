import matplotlib.pyplot as plt
import numpy as np

# Aircraft parameters

b_w = 20.4
L = b_w / 2

# Moment diagram

z_vals = np.arange(0, L, 0.1)
y_vals = np.sin(z_vals)
line = np.zeros_like(z_vals)

# Graph
plt.plot(z_vals, y_vals, label='Moment', color='red')
plt.plot(z_vals, line, linestyle='-', color='black')
plt.title('Moment Diagram')
plt.xlabel('z [m]')
plt.ylabel('M [Nm]')
plt.grid(True)
plt.xlim(0,L)

plt.legend()
plt.show()



