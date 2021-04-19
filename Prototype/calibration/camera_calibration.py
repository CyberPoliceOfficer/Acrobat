import numpy as np
import cv2 as cv
import glob

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

calib_grid_x = 9
calib_grid_y = 6
grid_side = 0.024 # grid square size in meter

# termination criteria
criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 30, 0.001)
# prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((calib_grid_x*calib_grid_y,3), np.float32)
objp[:,:2] = grid_side*np.mgrid[0:calib_grid_x,0:calib_grid_y].T.reshape(-1,2)

# Arrays to store object points and image points from all the images.
objpoints = [] # 3d point in real world space
imgpoints = [] # 2d points in image plane.
images = glob.glob('calibration_data/*.jpg')
for fname in images:
    img = cv.imread(fname)
    gray = cv.cvtColor(img, cv.COLOR_BGR2GRAY)
    # Find the chess board corners
    ret, corners = cv.findChessboardCorners(gray, (calib_grid_x,calib_grid_y), None)
    # If found, add object points, image points (after refining them)
    if ret == True:
        objpoints.append(objp)
        corners2 = cv.cornerSubPix(gray,corners, (11,11), (-1,-1), criteria)
        imgpoints.append(corners)
        # Draw and display the corners
        cv.drawChessboardCorners(img, (calib_grid_x,calib_grid_y), corners2, ret)
        
        cv.imshow('img',ResizeWithAspectRatio(img, width=1100))
        cv.waitKey(500) 
cv.destroyAllWindows()

# Calibrate the camera
ret, mtx, dist, rvecs, tvecs = cv.calibrateCamera(objpoints, imgpoints, gray.shape[::-1], None, None)

# Compute the mean error
mean_error = 0
for i in range(len(objpoints)):
    imgpoints2, _ = cv.projectPoints(objpoints[i], rvecs[i], tvecs[i], mtx, dist)
    error = cv.norm(imgpoints[i], imgpoints2, cv.NORM_L2)/len(imgpoints2)
    mean_error += error
print( "total error: {}".format(mean_error/len(objpoints)) )

# Save the calibration parameters
np.savez('cali_values', ret=ret, mtx=mtx, dist=dist, rvecs=rvecs, tvecs=tvecs)
