import numpy as np
import cv2 as cv
import glob

# [TESTING] Detect aruco markers

# Load the marker images
images = glob.glob('markers_data/*.jpg')
for fname in images:
   
    img = cv.imread(fname)
    aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
    aruco_side = 0.026

    cv.imshow('img', ResizeWithAspectRatio(img, width=1100))
    cv.waitKey()

    # Undistort the image
    img = cv.undistort(img, mtx, dist)
    img_gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)

    # Detect aruco markers
    corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img_gray, aruco_dict)
    rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)

    print(_objPoints)
    print(rvecs)
    print(tvecs)

    if (len(ids) > 0):
        cv.aruco.drawDetectedMarkers(img, corners, ids)
        rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
        for i in range(len(ids)):
            img = cv.aruco.drawAxis(img, mtx, dist, rvecs[i], tvecs[i],  0.05)
    cv.imshow('img', ResizeWithAspectRatio(img, width=1100))
    cv.waitKey()