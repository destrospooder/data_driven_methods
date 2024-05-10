#!/bin/env python3

'''
The Cooper Union
ME493 - Methods of Data-Driven Control
Prof. Dirk M. Luchtenburg

Final - MPC Implementation on a Frankenstein Predator-Prey Model
I guess this just commits genocide?
Benjamin Aziel
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import do_mpc 
from lotka_volterra import *

mpl.rcParams['font.size'] = 18
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['axes.grid'] = True

model_type = 'continuous'
model = do_mpc.model.Model(model_type)

x = model.set_variable(var_type = '_x', var_name = 'x', shape = (1, 1))
y = model.set_variable(var_type = '_x', var_name = 'y', shape = (1, 1))
u = model.set_variable(var_type = '_u', var_name = 'u', shape = (1, 1))

model.set_rhs('x', coeffs[0, :5] @ [x, y, u, x**2, x*y])
model.set_rhs('y', coeffs[1, :5] @ [x, y, u, x**2, x*y])

model.setup()

mpc = do_mpc.controller.MPC(model)

setup_mpc = {
    'n_horizon': 20,
    't_step': 0.1,
    'n_robust': 1,
    'store_full_solution': True,
}
mpc.set_param(**setup_mpc)

meyer = x**2 + y**2
lagrange = x**2 + y**2
mpc.set_objective(mterm = meyer, lterm = lagrange)

mpc.setup()

# controller sim
simulator = do_mpc.simulator.Simulator(model)
simulator.set_param(t_step = 0.1)
simulator.setup()

x0 = np.array([[100], [100]])
simulator.x0 = x0
mpc.x0 = x0

mpc.set_initial_guess()
mpc_graphics = do_mpc.graphics.Graphics(mpc.data)
sim_graphics = do_mpc.graphics.Graphics(simulator.data)

fig, ax = plt.subplots(2, figsize=(16,9))
fig.align_ylabels()

# capture
for g in [sim_graphics, mpc_graphics]:
    # Plot the angle positions (phi_1, phi_2, phi_2) on the first axis:
    g.add_line(var_type='_x', var_name='x', axis=ax[0])
    g.add_line(var_type='_x', var_name='y', axis=ax[0])

    # Plot the set motor positions (phi_m_1_set, phi_m_2_set) on the second axis:
    g.add_line(var_type='_u', var_name='u', axis=ax[1])

ax[0].set_ylabel('population')
ax[1].set_ylabel('unspecified control action')
ax[1].set_xlabel('time [yr]')

u0 = np.array([[0]])
for i in range(20):
    simulator.make_step(u0)

sim_graphics.plot_results()
# Reset the limits on all axes in graphic to show the data.
sim_graphics.reset_axes()

simulator.reset_history()
simulator.x0 = x0
mpc.reset_history()

for i in range(20):
    u0 = mpc.make_step(x0)
    x0 = simulator.make_step(u0)

# Plot predictions from t=0
mpc_graphics.plot_predictions(t_ind=0)
# Plot results until current time
sim_graphics.plot_results()
sim_graphics.reset_axes()
plt.legend()
plt.show()

print('debug')