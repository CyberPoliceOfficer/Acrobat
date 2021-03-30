from pololu_serial_master import servo_serial_master
from RSS_cc import RSS_cc
import numpy as np
import time

#serial_interface = servo_serial_master()

Manipulator = RSS_cc(r_b = 0.05257, r_m = 0.04814, d_b = 0.017, d_m = 0.08, d = 0.1175, h = 0.027)
alpha = Manipulator.inverse_kinematics([0, 0, 0.12, 0, 0, 0])
print(Manipulator.inverse_kinematic_jacobian([0, 0.0, 0.12, 0, 0, 0], alpha))

'''
alphas = Manipulator.inverse_kinematics([0, -0.04, 0.12, 0, 0, 0])
serial_interface.set_target_rad(alphas)
time.sleep(1.7)

alphas = Manipulator.inverse_kinematics([0, 0, 0.12, 0, 0, 0])
serial_interface.set_target_rad(alphas)
time.sleep(1.7)

alphas = Manipulator.inverse_kinematics([0.04, 0, 0.12, 0, 0, 0])
serial_interface.set_target_rad(alphas)
time.sleep(1.7)

alphas = Manipulator.inverse_kinematics([-0.04, 0, 0.12, 0, 0, 0])
serial_interface.set_target_rad(alphas)
time.sleep(1.7)

alphas = Manipulator.inverse_kinematics([0, 0, 0.12, 0, 0, 0])
serial_interface.set_target_rad(alphas)
time.sleep(1.7)

#serial_interface.set_target_rad(0*np.ones(6))
'''





'''
serial_interface.set_target_rad(2.9*np.ones(6))
time.sleep(0.7)
serial_interface.set_target_rad(-2.9*np.ones(6))
'''

'''
arr = np.array([0.7, -0.7, 0.7, -0.7, 0.7, -0.7])
serial_interface.set_target_rad(2*arr)
time.sleep(0.7)
serial_interface.set_target_rad(-2*arr)
time.sleep(0.7)
serial_interface.set_target_rad(2*arr)
time.sleep(0.7)
serial_interface.set_target_rad(-2*arr)
time.sleep(0.7)
serial_interface.set_target_rad(2*arr)
time.sleep(0.7)
serial_interface.set_target_rad(np.zeros(6))            
#serial_interface.set_target_rad(2.9*np.ones(6))
'''