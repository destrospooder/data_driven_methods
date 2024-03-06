#!/bin/env python

# just moving all my imported variables from the airfoil dataset here 
# so it's easier to grab em for dmd/sindy later

import h5py
import numpy as np

grid_file = h5py.File('dataset/airfoilDNS_grid.h5', 'r+')
data_file = h5py.File('dataset/airfoilDNS_a25f0p05.h5', 'r+')

# grid coords
x = np.squeeze(grid_file['x'][()])
y = np.squeeze(grid_file['y'][()])
nx = len(x)
ny = len(y)

# timesteps/timespan
t_field = np.squeeze(data_file['t_field'][()])
t_force = np.squeeze(data_file['t_force'][()])
nt = len(t_field)

# streamwise/transverse velocity
ux = np.squeeze(data_file['ux'][()])
uy = np.squeeze(data_file['uy'][()])

# airfoil coords
xa = np.squeeze(data_file['xa'][()]) 
ya = np.squeeze(data_file['ya'][()])