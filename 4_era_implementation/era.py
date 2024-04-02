#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import control as ct
import scipy as sc

k = 5
omega = np.array([2.4, 2.6, 6.5, 8.3, 9.3]) * 10**3 * 2 * np.pi
zeta = [0.03, 0.03, 0.042, 0.025, 0.032]
tau = 10e-4

s = ct.tf('s')

G_num_1 = k * omega[1]**2 * omega[2]**2 * omega[4]**2
G_num_2 = s**2 + 2 * zeta[0] * omega[0] * s + omega[0]**2
G_num_3 = s**2 + 2 * zeta[3] * omega[3] * s + omega[3]**2

G_den_1 = omega[0]**2 * omega[3]**2
G_den_2 = s**2 + 2 * zeta[1] * omega[1] * s + omega[1]**2
G_den_3 = s**2 + 2 * zeta[2] * omega[2] * s + omega[2]**2
G_den_4 = s**2 + 2 * zeta[4] * omega[4] * s + omega[4]**2

# the ctrl systems library doesn't support time delays directly, but we can use a pade approximation!
# (hell of a time to learn these existed)

delay = ct.pade(tau, 3)
G = G_num_1 * G_num_2 * G_num_3 * ct.tf(delay[0], delay[1]) / (G_den_1 * G_den_2 * G_den_3 * G_den_4)

# G = 1 / (s + 1)

T = np.linspace(0, 0.015, 1000)
t, y = ct.impulse_response(G, T)

plt.plot(t, y)
plt.xlabel('Time')
plt.ylabel('Output')
plt.title('Impulse Response for Atomic Force Microscope TF')
plt.tight_layout()
plt.savefig("era_figs/original_sys_impulse.png")
plt.show()

mc = int(np.floor((y.shape[0] - 1) / 2))
mo = mc

hank_schrader = np.zeros((mc, mo))
marie_schrader = np.zeros((mc, mo))

row_indices, col_indices = np.indices((mc, mo))
hank_schrader = y[row_indices + col_indices]
marie_schrader = y[row_indices + col_indices + 1]

# storage is cheaper than compute :)
try:
    U = np.load("svd_store/svdeez_u.npy")
    S = np.load("svd_store/svdeez_s.npy")
    Vh = np.load("svd_store/svdeez_vh.npy")
except FileNotFoundError:
    print("Performing SVD...")
    U, S, Vh = np.linalg.svd(hank_schrader, full_matrices = False)
    np.save("svd_store/svdeez_u.npy", U)
    np.save("svd_store/svdeez_s.npy", S)
    np.save("svd_store/svdeez_vh.npy", Vh)
    print("SVD components stored for future use.")

r = 9 # np.shape(U)[0]

Ur = U[:, :r]
SVr = np.diag(S[:r])
Vr = (Vh.T)[:, :r]

A = sc.linalg.fractional_matrix_power(SVr, -0.5) @ Ur.T @ marie_schrader @ Vr @ sc.linalg.fractional_matrix_power(SVr, -0.5)
B = (sc.linalg.fractional_matrix_power(SVr, 0.5) @ (Vr.T)[:, 0]).reshape(-1, 1)
C = Ur[0, :] @ sc.linalg.fractional_matrix_power(SVr, 0.5)

plt.stem(np.arange(r), np.square(S[:r]), "o")
plt.xlabel("Index")
plt.ylabel("Squared Hankel Singular Values")
plt.suptitle("Significant Hankel SV's (Squared)", fontsize=12)
plt.semilogy()
plt.savefig("era_figs/hankel_squared_svs.png")
plt.show()

G_new = ct.ss(A, B, C, 0, T[1])
ct.bode_plot([G, G_new * T[1]], wrap_phase = True)
plt.legend(['Original TF', 'ERA Reconstruction'])
plt.suptitle('Frequency Response Comparison')
plt.tight_layout()
plt.savefig("era_figs/freq_responses_compared.png")
plt.show()
