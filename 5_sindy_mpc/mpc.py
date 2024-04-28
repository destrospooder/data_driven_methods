
'''
The Cooper Union
ME493 - Methods of Data-Driven Control
Prof. Dirk M. Luchtenburg

Doing MPC on SINDyC-Derived System
Benjamin Aziel
'''

import numpy as np
import do_mpc as mpc
from sindyc import coeffs

model = mpc.model.Model('continuous')

x = model.set_variable(var_type = '_x', var_name = 'x', shape = (1, 1))
v = model.set_variable(var_type = '_x', var_name = 'v', shape = (1, 1))
F = model.set_variable(var_type = '_u', var_name = 'F', shape = (1, 1))

a11 = coeffs[0, 0]
a12 = coeffs[0, 1]
a21 = coeffs[1, 0]
a22 = coeffs[1, 1]

b1 = coeffs[0, 2]
b2 = coeffs[1, 2]

model.set_rhs('x', a11 * x + a12 * v + b1 * F)
model.set_rhs('v', a21 * x + a22 * v + b2 * F)

model.setup()

print('debug')