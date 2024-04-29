
'''
The Cooper Union
ME493 - Methods of Data-Driven Control
Prof. Dirk M. Luchtenburg

Implementing LQR using SINDyC-Derived Model
Benjamin Aziel
'''

import numpy as np
import matplotlib.pyplot as plt
import control as ctrl
from sindyc import coeffs, m, b, k, data, noise, t_eval

A = np.array([[0, 1], [-k/m, -b/m]])
B = np.array([[0], [1/m]])

A_sindy = coeffs[:2, :2]
B_sindy = np.vstack([coeffs[0, 2], coeffs[1, 2]])
print(B_sindy)
Q = np.eye(2)
R = 1

K, X, E = ctrl.lqr(A_sindy, B_sindy, Q, R)

H = ctrl.ss(A - B @ K, B @ K, np.eye(2), 0)

ctrl_t, ctrl_y = ctrl.initial_response(H, T = t_eval, X0 = [1, 0])

fig, ax = plt.subplots(1, 2)

ax[0].scatter(data.t, data.y[0] + noise[0, :], s=2)
ax[0].plot(ctrl_t, ctrl_y[0, :], color = 'red')
ax[0].set_title('Position vs. Time (Gaussian Noise, StdDev 1e-2) w/ Ctrl')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Position [m]')

ax[1].scatter(data.t, data.y[1] + noise[1, :], s=2)
ax[1].plot(ctrl_t, ctrl_y[1, :], color = 'red')
ax[1].set_title('Velocity vs. Time (Gaussian Noise, StdDev 1e-2) w/ Ctrl')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Velocity [m/s]')

plt.show()

print('debug')