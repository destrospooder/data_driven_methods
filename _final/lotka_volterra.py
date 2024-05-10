#!/bin/env python3

'''
The Cooper Union
ME493 - Methods of Data-Driven Control
Prof. Dirk M. Luchtenburg

Final - SINDYc Implementation on a Frankenstein Predator-Prey Model
Benjamin Aziel
'''

import numpy as np
import matplotlib.pyplot as plt
import pysindy as ps
from scipy.integrate import solve_ivp

noise_toggle = False
plots_toggle = False

# Lotka-Volterra predator-prey model
alpha = 1.1 # max prey per capita growth rate
beta = 0.4 # effect of the presence of predators on the prey death rate
gamma = 0.1 # predator per capita death rate
delta = 0.4 # effect of the presence of prey on predator birth rate

def u(t): # artificial inflation of the prey, let's say
    return np.sin(t)

def lotka_volterra(t, x, u):
    a = 0.1
    F = a * u(t)
    return(np.array([alpha * x[0] - beta * x[0] * x[1],
                    delta * x[0] * x[1] - gamma * x[1] + F]))

x0 = [10, 10]
t_span = [0, 80]
t_eval = np.linspace(t_span[0], t_span[1], 10000)

data = solve_ivp(lotka_volterra, 
                 t_span = t_span, 
                 y0 = x0, 
                 t_eval = t_eval,
                 args = (u,))

ctrl_inputs = u(t_eval)

fig, ax = plt.subplots(1, 2)

ax[0].scatter(data.t, data.y[0])
ax[0].set_title('Pop Density of Prey vs. Time (Noiseless)')
ax[0].set_xlabel('Time [yr]')
ax[0].set_ylabel('Population')

ax[1].scatter(data.t, data.y[1])
ax[1].set_title('Pop Density of Predators vs. Time (Noiseless)')
ax[1].set_xlabel('Time [yr]')
ax[1].set_ylabel('Population')

if plots_toggle: plt.show()

# just for funsies
noise_toggle = False

if noise_toggle:
    rng = np.random.seed(488) 
    noise = np.random.normal(0, 5e-2, np.shape(data.y))

    fig, ax = plt.subplots(1, 2)

    ax[0].scatter(data.t, data.y[0] + noise[0, :])
    ax[0].set_title('Pop Density of Prey vs. Time (Gaussian Noise, StdDev 1e-2)')
    ax[0].set_xlabel('Time [yr]')
    ax[0].set_ylabel('Population')

    ax[1].scatter(data.t, data.y[1] + noise[1, :])
    ax[1].set_title('Pop Density of Predators vs. Time (Gaussian Noise, StdDev 1e-2)')
    ax[1].set_xlabel('Time [yr]')
    ax[1].set_ylabel('Population')

if plots_toggle: plt.show()

else:
    noise = np.zeros(np.shape(data.y))

knob = 0.05

differentiation_method = ps.FiniteDifference(order = 2)
feature_library = ps.PolynomialLibrary(degree = 2, include_bias = False)
optimizer = ps.optimizers.STLSQ(threshold = knob)

model = ps.SINDy(
    differentiation_method = differentiation_method,
    feature_library = feature_library,
    optimizer = optimizer,
    feature_names = ['x', 'y', 'u']
)

X = np.vstack([data.y[0] + noise[0, :], data.y[1] + noise[1, :]])
model.fit(X.T, u = ctrl_inputs, t = data.t)
model.print()

coeffs = model.coefficients()
