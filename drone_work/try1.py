import time
import numpy as np

import cflib.crtp
from cflib.crazyflie.swarm import CachedCfFactory
from cflib.crazyflie.swarm import Swarm
from cflib.crazyflie.syncCrazyflie import SyncCrazyflie


def activate_led_bit_mask(scf):
    scf.cf.param.set_value('led.bitmask', 255)

def deactivate_led_bit_mask(scf):
    scf.cf.param.set_value('led.bitmask', 0)

def light_check(scf):
    activate_led_bit_mask(scf)
    time.sleep(2)
    deactivate_led_bit_mask(scf)
    time.sleep(2)

def take_off(scf):
    commander= scf.cf.high_level_commander

    commander.takeoff(1.0, 2)
    time.sleep(3)

def land(scf):
    commander= scf.cf.high_level_commander

    commander.land(0.0, 2)
    time.sleep(2)

    commander.stop()

def hover_sequence(scf):
    take_off(scf)
    land(scf)

def go_to_z(scf, x, y):
    commander = scf.cf.high_level_commander
    commander.go_to(x, y, 1, 0, 2)

def go_to_z_wrapper(scf, k, x_coords, y_coords):
    go_to_z(scf, x_coords[k], y_coords[k])

def run_sequence(scf: SyncCrazyflie, sequence):
    cf = scf.cf

    for arguments in sequence:
        commander = scf.cf.high_level_commander

        x, y, z = arguments[0], arguments[1], arguments[2]
        duration = arguments[3]

        print('Setting position {} to cf {}'.format((x, y, z), cf.link_uri))
        commander.go_to(x, y, z, 0, duration, relative=False)
        time.sleep(duration)

def dict_to_sorted_list(d):
    sorted_keys = sorted(d.keys())
    return [(key, d[key]) for key in sorted_keys]

uris = [
    'radio://0/70/2M/E7E7E7E7E7',
    'radio://0/80/2M/E7E7E7E7E7',
    'radio://0/80/2M/E7E7E7E7F0',
    'radio://0/80/2M/E7E7E7E7F2',
    # Add more URIs if you want more copters in the swarm
]

n = len(uris)

if __name__ == '__main__':
    cflib.crtp.init_drivers()
    factory = CachedCfFactory(rw_cache='./cache')
    with Swarm(uris, factory=factory) as swarm:
        print('Connected to Crazyflies')
        swarm.parallel_safe(light_check)
        swarm.reset_estimators()
        initial_positions = swarm.get_estimated_positions()
        sorted_positions = dict_to_sorted_list(initial_positions)
        print("Sorted positions:", sorted_positions)
        
        positions_array = np.array([pos for _, pos in sorted_positions])
        x = positions_array[:, 0]
        y = positions_array[:, 1]

        # Calculate line of best fit
        x_bar, y_bar = np.mean(x), np.mean(y)
        num = np.sum((x - x_bar) * (y - y_bar))
        den = np.sum((x - x_bar)**2)
        slope = num / den
        intercept = y_bar - slope * x_bar

        # Calculate closest points on the line of best fit
        points_inter_x = (slope / (slope**2 + 1)) * (x / slope + y - intercept)
        points_inter_y = slope * points_inter_x + intercept

        print("Original positions:", list(zip(x, y)))
        print("Closest points on line of best fit:", list(zip(points_inter_x, points_inter_y)))

        # Take off
        swarm.parallel_safe(take_off)
        print(swarm.get_estimated_positions())

        # Move to best fit points
        sequences = [[(points_inter_x[i], points_inter_y[i], 1, 3.0)] for i in range(n)]
        seq_args = {uri: [sequence] for uri, sequence in zip(uris, sequences)}

        swarm.parallel_safe(run_sequence, args_dict=seq_args)
        
        # Land
        swarm.parallel_safe(land)