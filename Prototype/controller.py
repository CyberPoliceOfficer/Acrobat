from pololu_serial_master import servo_serial_master
from RSS_cc import RSS_cc
import numpy as np
import time
from Acquisition import flir_cameras
import _thread
import cv2 as cv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 

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
		
def path():
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
	
def go_to ():
	# Setup
	aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
	aruco_side = 0.026
	npzfile = np.load('calibration/cali_values.npz')
	mtx = npzfile['mtx']
	dist = npzfile['dist']
	ret = npzfile['ret']
	rvecs = npzfile['rvecs']
	tvecs = npzfile['tvecs']

	platform_markers = np.array([13,14,18,19], dtype=np.int32)
	
	x_axis = np.array([26,27,28,29], dtype=np.int32)
	y_axis = np.array([26,21,16,11], dtype=np.int32)
	Manipulator = RSS_cc(r_b = 0.05257, r_m = 0.04814, d_b = 0.1624, d_m = 0.0832, d = 0.1175, h = 0.027)
	#Populu_maestro = servo_serial_master(port = '/dev/ttyACM0', baud_rate = 115200)

	# Goto
	#pose = Manipulator.home_pose()
	#thetas = Manipulator.inverse_kinematics (pose)
	#Populu_maestro.set_target_rad (thetas)
	#time.sleep(1)

	# Acquistion
	# camera = flir_cameras()
	#img = camera.Acquire()
	#np.savez('test_img', img=img)
	npzfile = np.load('test_img.npz')
	img = npzfile['img']
	corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img, aruco_dict, mtx, dist)
	cv.aruco.drawDetectedMarkers(img, corners, ids)
	rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
	platform_idexes = np.in1d(ids, platform_markers)
	# Compute

	# Find axis
	x_base = np.squeeze(tvecs[np.in1d(ids, x_axis[-1])] - tvecs[np.in1d(ids, x_axis[0])])
	x_base = x_base/np.linalg.norm(x_base)
	y_base = np.squeeze(tvecs[np.in1d(ids, y_axis[-1])] - tvecs[np.in1d(ids, y_axis[0])])
	y_base = y_base/np.linalg.norm(y_base)
	z_base = np.cross(x_base, y_base)
	print("Test (should be ~0):", np.matmul(y_base,x_base.T))

	O = np.array([x_base,y_base,z_base]).T    
	
	x_platform = np.squeeze(tvecs[np.in1d(ids, 13)] - tvecs[np.in1d(ids, 19)])
	x_platform = x_platform/np.linalg.norm(x_platform)
	y_platform = np.squeeze(tvecs[np.in1d(ids, 14)] - tvecs[np.in1d(ids, 19)])
	y_platform = y_platform/np.linalg.norm(y_platform)
	z_platform = np.cross(x_platform, y_platform)
	print("Test (should be ~0):", np.matmul(y_platform,x_platform.T))

	T = O @ np.squeeze((np.average(tvecs[platform_idexes == True], axis=0)))
	B = O @ np.squeeze((np.average(tvecs[platform_idexes == False], axis=0)))
	T = T - B
	B = B - B

	if True:
		x_base = O @ x_base
		y_base = O @ y_base
		z_base = O @ z_base

		x_platform = O @ x_platform
		y_platform = O @ y_platform
		z_platform = O @ z_platform 

	print(T)

	ax = plt.figure().add_subplot(projection='3d')
	
	if True:
		ax.quiver(B[0], B[1], B[2], x_base[0], x_base[1], x_base[2], color='red')
		ax.quiver(B[0], B[1], B[2], y_base[0], y_base[1], y_base[2], color='green')
		ax.quiver(B[0], B[1], B[2], z_base[0], z_base[1], z_base[2], color='blue')
	

	if True:
		ax.quiver(T[0], T[1], T[2], x_platform[0], x_platform[1], x_platform[2], color='red')
		ax.quiver(T[0], T[1], T[2], y_platform[0], y_platform[1], y_platform[2], color='green')
		ax.quiver(T[0], T[1], T[2], z_platform[0], z_platform[1], z_platform[2], color='blue')
	
	ax.quiver(0,0,0,T[0], T[1], T[2],  color='yellow')

	plt.show()


	M = np.array([x_platform,y_platform,z_platform]).T

	# might be wrong, need to doublecheck later
	R = (M @ O.T)
	print(R)

	# Display
	for i in range(len(ids)):
		img = cv.aruco.drawAxis(img, mtx, dist, rvecs[i], tvecs[i],  0.05)
	cv.imshow('img', img)
	#camera.Close()
	cv.waitKey()

if __name__ == "__main__":
    #path()
	go_to()

