import numpy as np
import time
from RSS_cc import RSS_cc
from Acquisition import flir_cameras
from pololu_serial_master import servo_serial_master
from random import uniform
import time
import cv2 as cv

def ResizeWithAspectRatio(image, width=None, height=None, inter=cv.INTER_AREA):
    dim = None
    (h, w) = image.shape[:2]

    if width is None and height is None:
        return image
    if width is None:
        r = height / float(h)
        dim = (int(w * r), height)
    else:
        r = width / float(w)
        dim = (width, int(h * r))

    return cv.resize(image, dim, interpolation=inter)


def main():
    aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
    aruco_side = 0.026
    # Load camera calibration
    npzfile = np.load('calibration/cali_values.npz')
    mtx = npzfile['mtx']
    dist = npzfile['dist']
    ret = npzfile['ret']
    rvecs = npzfile['rvecs']
    tvecs = npzfile['tvecs']

    camera = flir_cameras()
    Manipulator = RSS_cc(r_b = 0.05257, r_m = 0.04814, d_b = 0.017, d_m = 0.008, d = 0.1175, h = 0.027)
    Populu_maestro = servo_serial_master(port = '/dev/ttyACM0', baud_rate = 115200)
    Populu_maestro.set_speed_rad_s(3*np.ones(6))
    imgs = []
    us = []
    num_points = 400
    corner = 0.03
    for k in range(num_points):
        isNaN = True
        while (isNaN):
            p = np.zeros(6)
            p[0] = uniform(-corner, corner)
            p[1] = uniform(-corner, corner)
            p[2] = uniform(-corner, corner)
            p[3] = uniform(-0.3, 0.3)
            p[4] = uniform(-0.3, 0.3)
            p[5] = uniform(-0.3, 0.3)
            p = p + Manipulator.home_pose()
            thetas = Manipulator.inverse_kinematics(p)
            isNaN = np.isnan(np.sum(thetas))

        u = Populu_maestro.set_target_rad(thetas)
        us.append(u)
        time.sleep(2)
        
        # goto
        print(k,p)
        img = camera.Acquire()
        imgs.append(img)
        
        '''
        img = cv.undistort(img, mtx, dist)
        corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img, aruco_dict, mtx, dist)
        #recovered_ids, other = cv.aruco.refineDetectedMarkers(img, aruco_dict, corners, ids, rejectedImgPoints, mtx, dist)
        #print(len(recovered_ids))
        rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)

        print(len(ids))
        
    
        if (len(ids) > 0):
            cv.aruco.drawDetectedMarkers(img, corners, ids)
            rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
            for i in range(len(ids)):
                img = cv.aruco.drawAxis(img, mtx, dist, rvecs[i], tvecs[i],  0.05)
        cv.imshow('img', ResizeWithAspectRatio(img, width=1100))
        cv.waitKey()
        '''
        

    np.savez('calibration_data', images=imgs, us=us)
    camera.Close()

if __name__ == '__main__':
    main()
