#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pysindy as ps
from var_import import *

# storage is cheaper than compute :)
try:
    U = np.load("svd_store/svdeez_u.npy")
    S = np.load("svd_store/svdeez_s.npy")
    Vh = np.load("svd_store/svdeez_vh.npy")
except FileNotFoundError:
    from pod_analysis import U, S, Vh

V = np.transpose(Vh)

n_tempamps = 6  # how many temporal amplitudes are we considering?
tempamps = []
for i in range(n_tempamps):
    tempamps.append(S[i] * V[:, i])

tempamps = np.transpose(tempamps)
tempamps[:, 4:] = -tempamps[
    :, 4:
]  # for matlab standardization, scaling of evects is irrelevant

differentiation_method = ps.FiniteDifference(order=2)
feature_library = ps.PolynomialLibrary(degree=2)
optimizer = ps.STLSQ(threshold=0.0135)

model = ps.SINDy(
    differentiation_method=differentiation_method,
    feature_library=feature_library,
    optimizer=optimizer,
    feature_names=["x1", "x2", "x3", "x4", "x5", "x6"],
)

model.fit(tempamps, t=t_field - 50)
model.print()

x0 = np.array(
    [
        tempamps[0, 0],
        tempamps[0, 1],
        tempamps[0, 2],
        tempamps[0, 3],
        tempamps[0, 4],
        tempamps[0, 5],
    ]
).T

t_test = np.linspace(0, 35, 1000)

sim = model.simulate(x0, t=t_test)

fig, ax1 = plt.subplots(2, 3, figsize=(12, 6))
for k in range(6):
    ax = ax1.flat[k]
    ax.plot(t_test, sim[:, k])

plt.show()
