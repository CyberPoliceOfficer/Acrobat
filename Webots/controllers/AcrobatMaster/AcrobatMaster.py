from controller import Robot, Motor
from RSS_cc import RSS_cc
import numpy as np

a = 0.181
K1 = 1e-03
K2 = 0
TIME_STEP = 1

theta = np.radians(np.array([0, 330, 120, 90, 240, 210]))
phi = np.radians(np.array([0, 60, 120, 180, 240, 300]))
gamma = np.radians(np.array([35.26439, 90, 35.26439, 90, 35.26439, 90]))
z = np.array([-1, 1, -1, 1, -1, 1])
w = np.array([1, 1, 1, 1, 1, 1])

FM = np.array([0, 0, 0.2, 0, 0, 0])

A = np.array([K1*np.cos(theta)*np.sin(gamma),
              K1*np.sin(theta)*np.sin(gamma),
              K1*np.cos(gamma),
              K1*(a/np.sqrt(3)*np.sin(phi)*np.cos(gamma) - a/np.sqrt(6)*z*np.sin(theta)*np.sin(gamma)) - K2*w*np.cos(theta)*np.sin(gamma),
              K1*(a/np.sqrt(6)*z*np.cos(theta)*np.sin(gamma) - a/np.sqrt(3)*np.cos(phi)*np.cos(gamma)) - K2*w*np.sin(theta)*np.sin(gamma),
              K1*(a/np.sqrt(3)*np.cos(phi)*np.sin(theta)*np.sin(gamma) - a/np.sqrt(3)*np.sin(phi)*np.cos(theta)*np.sin(gamma)) - K2*w*np.cos(gamma)]);


u = np.linalg.solve(A, FM)
x = np.sign(u)
ang_vel = np.sqrt(np.abs(u))
ang_vel = np.multiply(ang_vel,x)

robot = Robot()

Manipulator = RSS_cc()

while robot.step(TIME_STEP) != -1:
    for i in range(1,7):
        propeller = robot.getDevice("prop" + str(i))
        propeller.setPosition(float('inf'))
        propeller.setVelocity(ang_vel[i-1])
        servo = robot.getDevice("servo" + str(i))
        alphas = Manipulator.inverse_kinematics([0.0,0.0,0.7,0,0,0])
        servo.setPosition(alphas[i-1])