import h5py
import numpy as np
import matplotlib.pyplot as plt

grid_file = h5py.File('airfoilDNS/airfoilDNS_grid.h5', 'r+')
data_file = h5py.File('airfoilDNS/airfoilDNS_a25f0p05.h5', 'r+')

mean_correction = True

base_angle = 25 # deg
freqs = 0.05 # rad/s

tstep = 1 # timestep

x = np.squeeze(grid_file['x'][()])
y = np.squeeze(grid_file['y'][()])
nx = len(x)
ny = len(y)

t_field = np.squeeze(data_file['t_field'][()])
t_force = np.squeeze(data_file['t_force'][()])
nt = len(t_field)

ux = np.squeeze(data_file['ux'][()]) # streamwise velocity
ux_reshaped = ux.reshape(nx * ny, nt, order = 'F').copy()

uy = np.squeeze(data_file['uy'][()]) # transverse velocity
uy_reshaped = uy.reshape(nx * ny, nt, order = 'F').copy()

X = np.vstack([ux_reshaped, uy_reshaped])

if mean_correction:
    X_mean= np.mean(X, axis = 1)
    X = X - np.reshape(X_mean, (len(X_mean), 1)) * np.ones([1, nt])

try:
    U = np.load('airfoilDNS/svdeez_u.npy')
    S = np.load('airfoilDNS/svdeez_s.npy')
    Vh = np.load('airfoilDNS/svdeez_vh.npy')
except FileNotFoundError:
    print("Performing SVD...")
    U, S, Vh = np.linalg.svd(X, full_matrices = False)
    np.save('airfoilDNS/svdeez_u.npy', U)
    np.save('airfoilDNS/svdeez_s.npy', S)
    np.save('airfoilDNS/svdeez_vh.npy', Vh)
    print("SVD components stored for future use.")

# Plotting squared singular values
"""plt.semilogy(np.arange(len(S)-1), np.square(S[:-1]), 'o')
plt.xlabel('Index')
plt.ylabel('Squared Singular Values')
plt.show()"""

# Spatial modes
midpoint = U.shape[0] // 2
U_ux = U[:midpoint, :]
U_uy = U[midpoint:, :]

"""fig, axes = plt.subplots(3, 2, figsize=(12, 8))

for i in range(6):
    row, col = i // 2, i % 2
    ax = axes[row, col]
    ax.plot(np.arange(len(U_ux[:,i])), U_ux[:,i])
    ax.set_xlabel('Index')
    ax.set_ylabel('Value')
    ax.set_title(f'U_ux[{i}]')

plt.tight_layout()
plt.show()

fig, axes = plt.subplots(3, 2, figsize=(12, 8))

for i in range(6):
    row, col = i // 2, i % 2
    ax = axes[row, col]
    ax.plot(np.arange(len(U_uy[:,i])), U_uy[:,i])
    ax.set_xlabel('Index')
    ax.set_ylabel('Value')
    ax.set_title(f'U_uy[{i}]')

plt.tight_layout()
plt.show()"""

MM = 0.01
v = np.arange(-1, 1.1, 0.1)

data = np.transpose(U_uy[:, 5].reshape(nx, ny, order = 'F'))
xv, yv = np.meshgrid(x, y)

plt.contourf(xv, yv, data, MM * v)
plt.colorbar()
plt.show()
