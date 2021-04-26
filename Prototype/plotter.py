import numpy as np
import cv2 as cv
import glob

aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
aruco_side = 0.026
npzfile = np.load('calibration/cali_values.npz')
mtx = npzfile['mtx']
dist = npzfile['dist']
ret = npzfile['ret']
rvecs = npzfile['rvecs']
tvecs = npzfile['tvecs']
parameters =  cv.aruco.DetectorParameters_create()
#parameters.minCornerDistanceRate = 0.001
#parameters.minMarkerPerimeterRate = 0.01
#parameters.maxMarkerPerimeterRate = 5
parameters.cornerRefinementMethod = cv.aruco.CORNER_REFINE_CONTOUR
#parameters.minDistanceToBorder = 1
#parameters.minMarkerDistanceRate = 0.1
#parameters.adaptiveThreshWinSizeMax = 35
parameters.polygonalApproxAccuracyRate = 0.03 

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


videos = glob.glob('calibration data/*')
for fname in videos:
    print(fname)
    npzfile = np.load(fname)
    images = npzfile['images']
    for img in images:
         # Undistort the image
        img = cv.undistort(img, mtx, dist)

        # Detect aruco markers
        corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img, aruco_dict, parameters=parameters)

        if (len(ids) > 0):
            cv.aruco.drawDetectedMarkers(img, corners, ids)
            rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
            for i in range(len(ids)):
                img = cv.aruco.drawAxis(img, mtx, dist, rvecs[i], tvecs[i],  0.05)
            cv.imshow('img', ResizeWithAspectRatio(img, width=1100))
            cv.waitKey()