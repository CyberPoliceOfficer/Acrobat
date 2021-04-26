from pololu_serial_master import servo_serial_master
from RSS_cc import RSS_cc
import numpy as np
import time
from Acquisition import flir_cameras
import _thread

global flag
flag = True
fps = 60


class controller:
	def __init__(self, manipulator = None, actuator_interface = None):
		# Interfaces
		self.manipulator = manipulator
		self.actuators = actuator_interface
	

		# Create the state struct
		self.state = {}
		self.state['p'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)
		self.state['v'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)
		self.state['theta'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)
		self.state['omega'] = np.array([0, 0, 0, 0, 0, 0], dtype=float)

	def go_home(self):
		home_pose = self.manipulator.home_pose()	
		self.actuators.set_speed_rad_s(np.array([0, 0, 0, 0, 0, 0], dtype=float))	
		self.actuators.set_target_rad(self.manipulator.inverse_kinematics(home_pose))
		self.state['p'] = home_pose

	def path_linear (self, target_point, v):
		self.T = 1e-3
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
		self.delta_p = 5e-3
		AB = target_point - self.state['p'][0:3]
		norm_AB = np.linalg.norm(AB)
		if (norm_AB < 1e-3):
			return
		AB_n = AB/norm_AB
		self.state['v'][0:3] = v * AB_n

		k = 1
		while (np.linalg.norm(target_point-self.state['p'][0:3]) >= self.delta_p):
			self.state['omega'] = self.manipulator.inverse_kinematic_jacobian(self.state['p']) @ self.state['v']
			self.actuators.set_speed_rad_s(np.abs(self.state['omega']))
			print(k,self.state['omega'])
			self.state['p'][0:3] += AB_n*self.delta_p
			print(k,self.state['p'])
			self.state['theta'] = self.manipulator.inverse_kinematics (self.state['p'])
			print(k,self.state['theta'])
			self.actuators.set_target_rad(self.state['theta'] )
			k+=1
			time.sleep(self.delta_p/v)

		self.state['omega'].fill(0)
		self.state['theta'].fill(0)

def camera_thread(corner):
	camera = flir_cameras()
	imgs = []
	while flag:
		imgs.append(camera.Acquire())
		time.sleep(1/fps)

	np.savez('path_04', images=imgs, fps=fps)
	camera.Close()
	print("Video saved and camera closed")

def path_thread(corner, speed, laps):
	global flag
	x = corner
	y = corner

	Manipulator = RSS_cc(r_b = 0.05257, r_m = 0.04814, d_b = 0.017, d_m = 0.008, d = 0.1175, h = 0.027)
	Populu_maestro = servo_serial_master(port = '/dev/ttyACM0', baud_rate = 115200)
	cnt = controller(manipulator = Manipulator, actuator_interface = Populu_maestro)
	
	cnt.go_home()
	cnt.path_linear(np.array([x, y, Manipulator.home_pose()[2]]),speed)
	
	time.sleep(1)
	for i in range(laps):
		cnt.path_linear(np.array([-x, y, Manipulator.home_pose()[2]]),speed)
		cnt.path_linear(np.array([-x, -y, Manipulator.home_pose()[2]]),speed)
		cnt.path_linear(np.array([x, -y, Manipulator.home_pose()[2]]),speed)
		cnt.path_linear(np.array([x, y, Manipulator.home_pose()[2]]),speed)
	print("Path ended, returning home")
	flag = False
	time.sleep(1)	
	cnt.go_home()
		
def main():
	corner = 0.04
	speed = 0.2
	laps = 10
	_thread.start_new_thread (path_thread, (corner,speed,laps,))
	camera = flir_cameras()
	imgs = []
	while flag:
		imgs.append(camera.Acquire())
		time.sleep(1/fps)

	file_name = 'path_videos/' + str(corner) + '_' + str(speed) + '_' + str(laps) + '_' + str(fps)
	np.savez(file_name, images=imgs, fps=fps, corner=corner, speed=speed, laps=laps)
	camera.Close()
	print("Video saved and camera closed")
	
	

if __name__ == "__main__":
    main()
	


