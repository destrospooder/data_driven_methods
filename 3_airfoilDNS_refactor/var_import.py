#!/bin/env python

# just moving all my imported variables from the airfoil dataset here 
# so it's easier to grab em for dmd/sindy later
# honestly might use a yaml file later haha

import h5py
import numpy as np

grid_file = h5py.File('dataset/airfoilDNS_grid.h5', 'r+')
data_file = h5py.File('dataset/airfoilDNS_a25f0p05.h5', 'r+')
params_file = h5py.File('dataset/airfoilDNS_parameters.h5', 'r+')

# spatial grid coords
x = np.squeeze(grid_file['x'][()])
y = np.squeeze(grid_file['y'][()])
nx = len(x)
ny = len(y)

# streamwise/transverse velocity
ux = np.squeeze(data_file['ux'][()])
uy = np.squeeze(data_file['uy'][()])

# airfoil coords
xa = np.squeeze(data_file['xa'][()]) 
ya = np.squeeze(data_file['ya'][()])
t_field = np.squeeze(data_file['t_field'][()])
dt_field = np.squeeze(params_file['dt_field'][()])
t_force = np.squeeze(data_file['t_force'][()])
nt = len(t_field)
