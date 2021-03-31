from pololu_serial_master import servo_serial_master
from RSS_cc import RSS_cc
import numpy as np
import time

class controller:
	def __init__(self, manipulator = None, actuator_interface = None, sampling_time = 1e-3):
		# Interfaces
		self.manipulator = manipulator
		self.actuators = actuator_interface
		
		# Signal processing
		self.T = sampling_time

		# Create the state struct
		self.state = {}
		self.state['p'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)
		self.state['v'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)
		self.state['theta'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)
		self.state['omega'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)

	def go_home(self):
		home_pose = self.manipulator.home_pose()		
		self.actuators.set_target_rad(self.manipulator.inverse_kinematics(home_pose))
		self.state['p'] = home_pose

	def path_linear (self, target_point, v):
		AB = target_point - self.state['p'][0:3]
		norm_AB = np.linalg.norm(AB)
		if (norm_AB < 1e-3):
			return

		self.state['v'][0:3] = v * AB/norm_AB
		k = 1
		while (np.linalg.norm(target_point-self.state['p'][0:3]) > v*self.T):
			self.state['omega'] = self.manipulator.inverse_kinematic_jacobian(self.state['p']) @ self.state['v']
			self.actuators.set_speed_rad_s(np.abs(self.state['omega']))
			print(k,self.state['omega'])
			self.state['p'][0:3] += self.T*self.state['v'][0:3]
			print(k,self.state['p'])
			self.state['theta'] = self.manipulator.inverse_kinematics (self.state['p'])
			print(k,self.state['theta'])
			self.actuators.set_target_rad(self.state['theta'] )
			k+=1
			time.sleep(self.T)

		self.state['omega'].fill(0)
		self.state['theta'].fill(0)
	
	def path_linear_v2 (self, target_point, v):
		AB = target_point - self.state['p'][0:3]
		norm_AB = np.linalg.norm(AB)
		if (norm_AB < 1e-3):
			return

		self.state['v'][0:3] = v * AB/norm_AB
		k = 1
		while (np.linalg.norm(target_point-self.state['p'][0:3]) > v*self.T):
			self.state['omega'] = self.manipulator.inverse_kinematic_jacobian(self.state['p']) @ self.state['v']
			self.actuators.set_speed_rad_s(np.abs(self.state['omega']))
			print(k,self.state['omega'])
			self.state['p'][0:3] += self.T*self.state['v'][0:3]
			print(k,self.state['p'])
			self.state['theta'] = self.manipulator.inverse_kinematics (self.state['p'])
			print(k,self.state['theta'])
			self.actuators.set_target_rad(self.state['theta'] )
			k+=1
			time.sleep(self.T)

		self.state['omega'].fill(0)
		self.state['theta'].fill(0)

			

def main():
	Manipulator = RSS_cc(r_b = 0.05257, r_m = 0.04814, d_b = 0.017, d_m = 0.08, d = 0.1175, h = 0.027)
	Populu_maestro = servo_serial_master(port = '/dev/ttyACM0', baud_rate = 115200)
	cnt = controller(manipulator = Manipulator, actuator_interface = Populu_maestro)
	cnt.go_home()
	#Populu_maestro.set_target_rad(0*np.pi/2*np.ones(6))
	corner=0.03
	speed=1	
	cnt.path_linear(np.array([corner, corner, Manipulator.home_pose()[2]]),speed)
	time.sleep(1)
	while True:
		cnt.path_linear(np.array([corner, corner, Manipulator.home_pose()[2]]),speed)
		cnt.path_linear(np.array([corner, -corner, Manipulator.home_pose()[2]]),speed)
		cnt.path_linear(np.array([-corner, -corner, Manipulator.home_pose()[2]]),speed)
		cnt.path_linear(np.array([-corner, corner, Manipulator.home_pose()[2]]),speed)
	cnt.go_home()


if __name__ == "__main__":
    main()
	


