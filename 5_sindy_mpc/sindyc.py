
'''
The Cooper Union
ME493 - Methods of Data-Driven Control
Prof. Dirk M. Luchtenburg

SINDYc Implementation on a Mass-Spring-Damper System
Benjamin Aziel
'''

import numpy as np
import matplotlib.pyplot as plt
import pysindy as ps
from scipy.integrate import solve_ivp

m = 1 # kg
b = 0.4 # kg/s
k = 1 # N/m

def u(t):
    return 2

def mss(t, x, u):
    F = u(t)
    return(np.array([x[1], -k/m * x[0] - b/m * x[1] + F/m]))

x0 = [1, 0]
t_span = [0, 10]
t_eval = np.linspace(t_span[0], t_span[1], 10000)

data = solve_ivp(mss, 
                 t_span = t_span, 
                 y0 = x0, 
                 t_eval = t_eval,
                 args = (u,))

ctrl_inputs = u(t_eval)

fig, ax = plt.subplots(1, 2)

ax[0].scatter(data.t, data.y[0])
ax[0].set_title('Position vs. Time (Noiseless)')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Position [m]')

ax[1].scatter(data.t, data.y[1])
ax[1].set_title('Velocity vs. Time (Noiseless)')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Velocity [m/s]')

# plt.show()

# just for funsies
noise = np.random.normal(0, 1e-3, np.shape(data.y))

fig, ax = plt.subplots(1, 2)

ax[0].scatter(data.t, data.y[0] + noise[0, :])
ax[0].set_title('Position vs. Time (Gaussian Noise, StdDev 1e-3)')
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel('Position [m]')

ax[1].scatter(data.t, data.y[1] + noise[1, :])
ax[1].set_title('Velocity vs. Time (Gaussian Noise, StdDev 1e-3)')
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel('Velocity [m/s]')

# plt.show()

knob = 0.75

differentiation_method = ps.FiniteDifference(order = 2)
feature_library = ps.PolynomialLibrary(degree = 1, include_bias = False)
optimizer = ps.optimizers.STLSQ(threshold = knob)

model = ps.SINDy(
    differentiation_method = differentiation_method,
    feature_library = feature_library,
    optimizer = optimizer,
    feature_names = ['x', 'v', 'F']
)

X = np.vstack([data.y[0] + noise[0, :], data.y[1] + noise[1, :]])
model.fit(X.T, u = ctrl_inputs, t = data.t)
model.print()

coeffs = model.coefficients()
