#!/bin/env python

import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from var_import import *

# Doing the manual implementation first using Brunton example code, then I'll learn how to use PySINDy


def poolData(yin, nVars, polyorder):
    n = yin.shape[0]
    yout = np.zeros((n, 1))

    # poly order 0
    yout[:, 0] = np.ones(n)

    # poly order 1
    for i in range(nVars):
        yout = np.append(yout, yin[:, i].reshape((yin.shape[0], 1)), axis=1)

    # poly order 2
    if polyorder >= 2:
        for i in range(nVars):
            for j in range(i, nVars):
                yout = np.append(
                    yout, (yin[:, i] * yin[:, j]).reshape((yin.shape[0], 1)), axis=1
                )

    # poly order 3
    if polyorder >= 3:
        for i in range(nVars):
            for j in range(i, nVars):
                for k in range(j, nVars):
                    yout = np.append(
                        yout,
                        (yin[:, i] * yin[:, j] * yin[:, k]).reshape((yin.shape[0], 1)),
                        axis=1,
                    )

    return yout


def sparsifyDynamics(Theta, Xdot, sparse_knob, n):
    # sequentially thresholded least-squares, see dynamicslab github
    Xi = np.linalg.lstsq(Theta, Xdot, rcond=None)[0]

    for k in range(10):
        smallinds = np.abs(Xi) < sparse_knob
        # where the sparse_knob is represented by lambda in databookV2

        Xi[smallinds] = 0

        for ind in range(n):
            biginds = smallinds[:, ind] == 0
            Xi[biginds, ind] = np.linalg.lstsq(
                Theta[:, biginds], Xdot[:, ind], rcond=None
            )[0]

    return Xi


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
dtempamps = np.diff(tempamps, axis=0) / dt_field
tempamps = tempamps[:-1]

Theta = poolData(tempamps, n_tempamps, 2)  # because fluids is a quadratic domain
sparse_knob = 0.01  # 0.01381477
Xi = sparsifyDynamics(Theta, dtempamps, sparse_knob, n_tempamps)
print(Xi)

fdot_sindy = Theta @ Xi

f_sindy = np.hstack(
    [cumtrapz(fdot_sindy[:, k], t_field[:-1])[:, None] for k in range(6)]
)

fig, ax1 = plt.subplots(2, 3, figsize=(12, 6))

for k in range(6):
    ax = ax1.flat[k]
    ax.plot(t_field[:-2], f_sindy[:, k] + tempamps[0, k])
    ax.plot(t_field[:-1], tempamps[:, k], "r-")

    fig.suptitle(
        f"sindy-derived temporal amplitudes (lambda = {float(sparse_knob)}) vs. pod temporal amplitudes",
        fontsize=12,
    )
    plt.legend(
        [f"SINDy (lambda = {float(sparse_knob)})", "POD"],
        bbox_to_anchor=(1.05, 0),
        loc="lower left",
    )
    plt.tight_layout()
    plt.savefig("sindy_figs/sindy_fit.png")

plt.show()
