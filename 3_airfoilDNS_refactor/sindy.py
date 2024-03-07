#!/bin/env python

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from var_import import dt_field, t_field

# Doing the manual implementation first using Brunton example code, then I'll learn how to use PySINDy

def poolData(yin,nVars,polyorder):
    n = yin.shape[0]
    yout = np.zeros((n,1))
    
    # poly order 0
    yout[:,0] = np.ones(n)
    
    # poly order 1
    for i in range(nVars):
        yout = np.append(yout,yin[:,i].reshape((yin.shape[0],1)),axis=1)
    
    # poly order 2
    if polyorder >= 2:
        for i in range(nVars):
            for j in range(i,nVars):
                yout = np.append(yout,(yin[:,i]*yin[:,j]).reshape((yin.shape[0],1)),axis=1)
                
    # poly order 3
    if polyorder >= 3:
        for i in range(nVars):
            for j in range(i,nVars):
                for k in range(j,nVars):
                    yout = np.append(yout,(yin[:,i]*yin[:,j]*yin[:,k]).reshape((yin.shape[0],1)),axis=1)
    
    return yout

def sparsifyDynamics(Theta, Xdot, sparse_knob, n):
# sequentially thresholded least-squares, see dynamicslab github
    Xi = np.linalg.lstsq(Theta, Xdot, rcond = None)[0]

    for k in range(10):
        smallinds = np.abs(Xi) < sparse_knob
        # where the sparse_knob is represented by lambda in databookV2

        Xi[smallinds] = 0

        for ind in range(n):
            biginds = smallinds[:, ind] == 0
            Xi[biginds, ind] = np.linalg.lstsq(Theta[:, biginds], Xdot[:, ind], rcond = None)[0]

    return Xi

# storage is cheaper than compute :)
try:
    U = np.load('svd_store/svdeez_u.npy')
    S = np.load('svd_store/svdeez_s.npy')
    Vh = np.load('svd_store/svdeez_vh.npy')
except FileNotFoundError:
    from pod_analysis import U, S, Vh

V = np.transpose(Vh)

n_tempamps = 6 # how many temporal amplitudes are we considering?
tempamps = []
for i in range(n_tempamps):
    tempamps.append(S[i] * V[:, i])

tempamps = np.transpose(tempamps)
dtempamps = np.diff(tempamps, axis = 0) / dt_field
dtempamps = np.vstack([dtempamps, dtempamps[399, :]])

Theta = poolData(tempamps, n_tempamps, 2) # because fluids is a quadratic domain
sparse_knob = 0.0135
Xi = sparsifyDynamics(Theta, dtempamps, sparse_knob, n_tempamps)
print(Xi)

def idd_system(t, x):
    # x is a state vector consisting of x1, x2, x3, ...
    theta = np.array([
             1, 
             x[0], 
             x[1], 
             x[2], 
             x[3], 
             x[4], 
             x[5], 
             x[0] * x[0],
             x[0] * x[1],
             x[0] * x[2],
             x[0] * x[3],
             x[0] * x[4],
             x[0] * x[5],
             x[1] * x[1],
             x[1] * x[2],
             x[1] * x[3],
             x[1] * x[4],
             x[1] * x[5],
             x[2] * x[2],
             x[2] * x[3],
             x[2] * x[4],
             x[2] * x[5],
             x[3] * x[3],
             x[3] * x[4],
             x[3] * x[5],
             x[4] * x[4],
             x[4] * x[5],
             x[5] * x[5]])

    return np.array([
        Xi[:, 0] @ theta,
        Xi[:, 1] @ theta,
        Xi[:, 2] @ theta,
        Xi[:, 3] @ theta,
        Xi[:, 4] @ theta,
        Xi[:, 5] @ theta,
    ])
 
x0 = np.array([tempamps[0, 0],
               tempamps[0, 1],
               tempamps[0, 2],
               tempamps[0, 3],
               tempamps[0, 4],
               tempamps[0, 5]])

t = [0, 40]
t_eval = np.linspace(t[0], t[1], 1000) # just to beef up resolution of solve_ivp

sol = solve_ivp(idd_system, t, x0, t_eval = t_eval)
print(sol)

fig, ax1 = plt.subplots(2, 3, figsize = (12, 6))

for k in range(6):
    ax = ax1.flat[k]
    ax.plot(sol.t + 50, sol.y[k, :])
    ax.plot(t_field, S[k] * V[:, k], 'r-')

fig.suptitle("sindy-derived temporal amplitudes (lambda = 0.0135) vs. pod temporal amplitudes", fontsize = 12)
plt.tight_layout()
plt.savefig('sindy_figs/comparison.png')
plt.show()