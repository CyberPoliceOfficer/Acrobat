import numpy as np
import cv2 as cv
import glob
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
aruco_side = 0.026
npzfile = np.load('calibration/cali_values.npz')
mtx = npzfile['mtx']
dist = npzfile['dist']
ret = npzfile['ret']
rvecs = npzfile['rvecs']
tvecs = npzfile['tvecs']
platform_markers = np.array([13,14,18,19], dtype=np.int32)

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

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

videos = glob.glob('path_videos/*')
for fname in videos:
    print(fname)
    npzfile = np.load(fname)
    images = npzfile['images']
    T = np.empty( shape=(0, 0) )
    for img in images:
        # Detect aruco markers
        corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img, aruco_dict, mtx, dist)
        
        platform_idexes = np.in1d(ids, platform_markers)
        if (np.sum(platform_idexes) == 4):
            rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
            if (T.size == 0):
                T = np.average(tvecs[platform_idexes == True], axis=0)
            else:
                T = np.vstack((T, np.average(tvecs[platform_idexes == True], axis=0)))

    print(np.shape(T))
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot (T[:,0], T[:,1], T[:,2])
    set_axes_equal(ax)
    

    plt.show()