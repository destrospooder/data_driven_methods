import h5py
import numpy as np
import matplotlib.pyplot as plt

grid_file = h5py.File('airfoilDNS_grid.h5', 'r+')
data_file = h5py.File('airfoilDNS_a25f0p05.h5', 'r+')

mean_correction = True

x = np.squeeze(grid_file['x'][()])
y = np.squeeze(grid_file['y'][()])
nx = len(x)
ny = len(y)

t_field = np.squeeze(data_file['t_field'][()])
t_force = np.squeeze(data_file['t_force'][()])
nt = len(t_field)

ux = np.squeeze(data_file['ux'][()]) # streamwise velocity
ux_reshaped = np.transpose(ux) # pains translating matlab to py code
ux_reshaped = ux_reshaped.reshape(nx * ny, nt, order = 'F')

uy = np.squeeze(data_file['uy'][()]) # transverse velocity
uy_reshaped = np.transpose(uy) # pains translating matlab to py code
uy_reshaped = uy_reshaped.reshape(nx * ny, nt, order = 'F')

xa = np.squeeze(data_file['xa'][()])
ya = np.squeeze(data_file['ya'][()])

X = np.vstack([ux_reshaped, uy_reshaped])

if mean_correction:
    X_mean= np.mean(X, axis = 1)
    X_mean_temp = np.transpose(X_mean) # pains translating matlab to py code
    X = X - np.reshape(X_mean_temp, (len(X_mean_temp), 1)) * np.ones([1, nt])

# storage is cheaper than compute
try:
    U = np.load('svd_store/svdeez_u.npy')
    S = np.load('svd_store/svdeez_s.npy')
    Vh = np.load('svd_store/svdeez_vh.npy')
except FileNotFoundError:
    print("Performing SVD...")
    U, S, Vh = np.linalg.svd(X, full_matrices = False)
    np.save('svd_store/svdeez_u.npy', U)
    np.save('svd_store/svdeez_s.npy', S)
    np.save('svd_store/svdeez_vh.npy', Vh)
    print("SVD components stored for future use.")

V = np.transpose(Vh)
SV = np.diag(S)

# Plotting squared singular values
plt.semilogy(np.arange(len(S) - 1), np.square(S[:-1]), 'o')
plt.xlabel('Index')
plt.ylabel('Squared Singular Values')
plt.suptitle("squared sv's", fontsize = 12)
plt.savefig('pod_analysis_figs/squared_sv.png')
plt.show()

# There's a pretty noticeable elbow around r = 16 so let's plot those in particular
plt.stem(np.arange(16), np.square(S[:16]), 'o')
plt.xlabel('Index')
plt.ylabel('Squared Singular Values')
plt.suptitle("first 16 squared sv's", fontsize = 12)
plt.semilogy()
plt.savefig('pod_analysis_figs/squared_sv_truncated.png')
plt.show()

U_ux = U[0:len(ux_reshaped), :]
U_uy = U[len(uy_reshaped):, :]

MM = 0.01
v = np.arange(-2, 2.1, 0.1)

# ux spatial modes
fig, ax1 = plt.subplots(2, 3, figsize = (12, 6))

for k in range(6):
    ax = ax1.flat[k]
    im = ax.contourf(x, y, np.transpose(U_ux[:, k].reshape(nx, ny, order = 'F')), levels = MM * v, vmin = -MM, vmax = MM)
    fig.colorbar(im, ax = ax)
    ax.set_title(f"ux mode {k + 1}")

fig.suptitle("ux spatial modes", fontsize = 12)
plt.tight_layout()
plt.savefig('pod_analysis_figs/ux_spatial_modes.png')
plt.show()

# uy spatial modes
fig, ax2 = plt.subplots(2, 3, figsize = (12, 6))

for k in range(6):
    ax = ax2.flat[k]
    im = ax.contourf(x, y, np.transpose(U_uy[:, k].reshape(nx, ny, order = 'F')), levels = MM * v, vmin = -MM, vmax = MM)
    fig.colorbar(im, ax = ax)
    ax.set_title(f"uy mode {k + 1}")
    # modes 5, 6 seem to have weird inversions in their behavior when compared to matlab. debug later

fig.suptitle("uy spatial modes", fontsize = 12)
plt.tight_layout()
plt.savefig('pod_analysis_figs/uy_spatial_modes.png')
plt.show()

# temporal amplitudes
fig, ax3 = plt.subplots(2, 3, figsize = (12, 6))

for k in range(6):
    ax = ax3.flat[k]
    ax.plot(t_field, S[k] * V[:, k])

fig.suptitle("temporal amplitudes", fontsize = 12)
plt.tight_layout()
plt.savefig('pod_analysis_figs/temporal_amplitudes.png')
plt.show()

# reconstruction for rank 4
r = 4
X_approx = U[:, :r] @ SV[:r, :r] @ Vh[:r, :]
t_star = np.linspace(0, nt/2 - 1, 6) # from temporal amplitudes we see that we get two periods out of this bad boy

# ux reconstruction
fig, ax4 = plt.subplots(2, 3, figsize = (12, 6))

cbar_max = np.array([])
cbar_min = np.array([])

for k in range(6):
    cbar_max = np.append(cbar_max, np.quantile(X_approx[0:nx*ny, int(t_star[k])], 0.99)) # outlier clipping
    cbar_min = np.append(cbar_min, np.quantile(X_approx[0:nx*ny, int(t_star[k])], 0.01))
cbar_max = np.max(1.01 * cbar_max)
cbar_min = np.min(1.01 * cbar_min)

for k in range(6):
    ax = ax4.flat[k]
    ux_approx = np.transpose(X_approx[0:nx*ny, int(t_star[k])]) # pains translating matlab to py code
    ux_approx = ux_approx.reshape(nx, ny, order = 'F')
    im = ax.contourf(x, y, np.transpose(ux_approx), levels = np.linspace(cbar_min, cbar_max, 40), vmin = cbar_min, vmax = cbar_max, extend = 'both')
    fig.colorbar(im, ax = ax)
    ax.set_title(f"timestep {int(t_star[k])}")
    ax.plot(xa, ya, 'r-') # plot all airfoil locations

fig.suptitle("ux reconstruction - rank 4", fontsize = 12)
plt.tight_layout()
plt.savefig('pod_analysis_figs/ux_reconstruction.png')
plt.show()

# uy reconstruction
fig, ax5 = plt.subplots(2, 3, figsize = (12, 6))

cbar_max = np.array([])
cbar_min = np.array([])

for k in range(6):
    cbar_max = np.append(cbar_max, np.quantile(X_approx[0:nx*ny, int(t_star[k])], 0.99)) # outlier clipping
    cbar_min = np.append(cbar_min, np.quantile(X_approx[0:nx*ny, int(t_star[k])], 0.01))

cbar_max = np.max(1.01 * cbar_max) 
cbar_min = np.min(1.01 * cbar_min)

for k in range(6):
    ax = ax5.flat[k]
    uy_approx = np.transpose(X_approx[nx*ny:, int(t_star[k])]) # pains translating matlab to py code
    uy_approx = uy_approx.reshape(nx, ny, order = 'F')
    im = ax.contourf(x, y, np.transpose(uy_approx), levels = np.linspace(cbar_min, cbar_max, 40), vmin = cbar_min, vmax = cbar_max, extend = 'both')
    fig.colorbar(im, ax = ax)
    ax.set_title(f"timestep {int(t_star[k])}")
    ax.plot(xa, ya, 'r-') # plot all airfoil locations

fig.suptitle("uy reconstruction - rank 4", fontsize = 12)
plt.tight_layout()
plt.savefig('pod_analysis_figs/uy_reconstruction.png')
plt.show()

# for debugging
print('heeheeheehaw')
print('hoohaahaahaa')