from controller import Robot, Motor, GPS
from RSS_cc import RSS_cc
import numpy as np
import scipy.io

TIME_STEP = 1

servos = []
robot = Robot()
print(robot.getCustomData)
gps = robot.getDevice("gps")
gps.enable(TIME_STEP)

Manipulator = RSS_cc()

time_stamp = np.empty((1,0))
poses = np.empty((3,0))

joint_speed = np.empty((6,0))
alpha_tmp = np.zeros((6,1))

while robot.step(TIME_STEP) != -1:
    t = robot.getTime()
    time_stamp = np.append(time_stamp, t)
    pose = np.array(gps.getValues()).reshape(3,1)
    poses = np.append(poses, pose, axis=1)
    for i in range(1,7):
        servo = robot.getDevice("servo" + str(i))
        alpha_tmp[i-1] = servo.getVelocity()
    
    joint_speed = np.append(joint_speed, alpha_tmp, axis=1)
    
    if (t<=1 and t>0):
        for i in range(1,7):
            servo = robot.getDevice("servo" + str(i))
            alphas = Manipulator.inverse_kinematics([0.1,0.1,0.7,0,0,0])
            servo.setPosition(alphas[i-1])
    if (t<=2 and t>1):
        for i in range(1,7):
            servo = robot.getDevice("servo" + str(i))
            alphas = Manipulator.inverse_kinematics([0.2,0.2,0.7,0,0,0])
            servo.setPosition(alphas[i-1])
    
    if (t > 20):
        break

           
            
    
         
#scipy.io.savemat('JacobianVal.mat', {'poses': poses, 'time_stamp': time_stamp, 'joint_speed': joint_speed})
print("Done!")

