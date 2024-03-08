#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from var_import import *

tv = 75 # truncation value

def calc_DMD(X, tv):
    X1 = X[:, :-1]
    X2 = X[:, 1:]

    U1, S1, Vh1 = np.linalg.svd(X1, full_matrices = False)

    U1 = U1[:, :tv]
    SV1 = np.diag(S1[:tv])
    Vh1 = Vh1[:tv, :]

    A_tilde = np.linalg.solve(SV1.T, (U1.T @ X2 @ Vh1.T).T).T
    Lambda, W = np.linalg.eig(A_tilde)
    Lambda = np.diag(Lambda)

    Phi = X2 @ np.linalg.solve(SV1.T, Vh1).T @ W
    b = np.linalg.solve(W @ Lambda, SV1 @ Vh1[:, 0])

    sort_inds = np.argsort(np.abs(b))[::-1]
    b = b[sort_inds]
    Lambda = np.diag(Lambda)
    Lambda = Lambda[sort_inds]
    Phi = Phi[:, sort_inds]

    return A_tilde, Phi, Lambda, b

# can probably put this in a class or something - that's a later problem though
mean_correction = True

ux_reshaped = np.transpose(ux) # pains translating matlab to py code
ux_reshaped = ux_reshaped.reshape(nx * ny, nt, order = 'F')

uy_reshaped = np.transpose(uy) # pains translating matlab to py code
uy_reshaped = uy_reshaped.reshape(nx * ny, nt, order = 'F')

X = np.vstack([ux_reshaped, uy_reshaped])

if mean_correction:
    X_mean= np.mean(X, axis = 1)
    X_mean_temp = np.transpose(X_mean) # pains translating matlab to py code
    X = X - np.reshape(X_mean_temp, (len(X_mean_temp), 1)) * np.ones([1, nt])

A_tilde, Phi, Lambda, b = calc_DMD(X, tv)

Eigscts = np.log(Lambda) / dt_field

plt.stem(np.imag(Eigscts) / (2 * np.pi), abs(b * Lambda) / np.max(abs(b * Lambda)), '-')
plt.xlim([-0.02, 1.5])
plt.ylim([0, 1.2])
plt.xlabel('frequency')
plt.ylabel('dmd mode amplitude (scaled)')
plt.suptitle(f'dmd mode amplitudes w/ truncation value {int(tv)}', fontsize = 12)
plt.savefig('dmd_analysis_figs/mode_amplitudes.png')
plt.show()

fig, ax1 = plt.subplots(2, 3, figsize = (12, 6))

MM = 0.01
v = np.arange(-1, 1.1, 0.1)

for k in range(6):
    ax = ax1.flat[k]
    k = 2 * k # remove cpx conjugate modes
    Phi_temp = Phi[int(len(Phi) / 2):, k].reshape(nx, ny, order = 'F')
    im = ax.contourf(x, y, np.transpose(Phi_temp), levels = MM * v, vmin = -MM, vmax = MM, extend = 'both')
    ax.plot(xa, ya, 'r-') # plot all airfoil locations
    ax.set_title(f'frequency = {float(abs(np.imag(Eigscts[k]) / (2 * np.pi))):.5f}')
    fig.colorbar(im, ax = ax)

fig.suptitle("uy dmd modes", fontsize = 12)
plt.tight_layout()
plt.savefig('dmd_analysis_figs/uy_dmd_modes.png')
plt.show()

fig, ax2 = plt.subplots(2, 3, figsize = (12, 6))

MM = 0.01
v = np.arange(-1, 1.1, 0.1)

for k in range(6):
    ax = ax2.flat[k]
    k = 2 * k # remove cpx conjugate modes
    Phi_temp = Phi[:int(len(Phi) / 2), k].reshape(nx, ny, order = 'F')
    im = ax.contourf(x, y, np.transpose(Phi_temp), levels = MM * v, vmin = -MM, vmax = MM, extend = 'both')
    ax.plot(xa, ya, 'r-') # plot all airfoil locations
    ax.set_title(f'frequency = {float(abs(np.imag(Eigscts[k]) / (2 * np.pi))):.5f}')
    fig.colorbar(im, ax = ax)

fig.suptitle("ux dmd modes", fontsize = 12)
plt.tight_layout()
plt.savefig('dmd_analysis_figs/ux_dmd_modes.png')
plt.show()
