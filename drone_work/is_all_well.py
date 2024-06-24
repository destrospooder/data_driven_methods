import logging
import sys
import time
from threading import Event

import cflib.crtp
from cflib.crazyflie import Crazyflie
from cflib.crazyflie.log import LogConfig
from cflib.crazyflie.syncCrazyflie import SyncCrazyflie
from cflib.positioning.position_hl_commander import PositionHlCommander
from cflib.utils import uri_helper

URI = uri_helper.uri_from_env(default='radio://0/80/2M/E7E7E7E7F0')

deck_attached_event = Event()

logging.basicConfig(level=logging.ERROR)
position_estimate = [0, 0, 0]

def log_callback(timestamp, data, logconf):
    # print(data)
    global position_estimate
    position_estimate[0] = data['stateEstimate.x']
    position_estimate[1] = data['stateEstimate.y']
    position_estimate[2] = data['stateEstimate.z']
    vbat = data['pm.vbat']
    battery_percentage = round(((vbat - 3.0) / (4.23 - 3.0)) * 100, 2)  # Calculate and round battery percentage
    print(f"x = {round(position_estimate[0], 3)} m, y = {round(position_estimate[1], 3)} m, z = {round(position_estimate[2], 3)} m")
    print(f"Battery voltage: {round(vbat, 3)} V, Battery percentage: {battery_percentage}%")

def param_deck_flow(_, value_str):
    value = int(value_str)
    if value:
        deck_attached_event.set()
        print('Deck is attached!')
    else:
        print('Deck is NOT attached!')

if __name__ == '__main__':
    cflib.crtp.init_drivers()

    with SyncCrazyflie(URI, cf=Crazyflie(rw_cache='./cache')) as scf:

        scf.cf.param.add_update_callback(group='deck', name='bcLighthouse4',
                                         cb=param_deck_flow)
        time.sleep(1)

        logconf = LogConfig(name='Status', period_in_ms=10)
        logconf.add_variable('stateEstimate.x', 'float')
        logconf.add_variable('stateEstimate.y', 'float')
        logconf.add_variable('stateEstimate.z', 'float')
        logconf.add_variable('pm.vbat', 'float')
        scf.cf.log.add_config(logconf)
        logconf.data_received_cb.add_callback(log_callback)

        if not deck_attached_event.wait(timeout=5):
                print('No flow deck detected!')
                sys.exit(1)

        logconf.start()
        time.sleep(1)
        logconf.stop()
