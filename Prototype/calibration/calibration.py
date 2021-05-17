import numpy as np
import cv2 as cv
import glob

def extract_points_from_video ():
    aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
    aruco_side = 0.026
    # Load camera calibration
    npzfile = np.load('cali_values.npz')
    mtx = npzfile['mtx']
    dist = npzfile['dist']
    ret = npzfile['ret']
    rvecs = npzfile['rvecs']
    tvecs = npzfile['tvecs']
    videos = glob.glob('RSS_calibration_data/*')
    ts = np.empty( shape=(0, 0) )
    platform_markers = np.array([13,14,18,19], dtype=np.int32)
    for fname in videos:
        print(fname)
        npzfile = np.load(fname)
        images = npzfile['images']
        us = npzfile['us']
        for img, us in zip(images, us):

            # Detect aruco markers
            corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img, aruco_dict, mtx, dist)
            
            platform_idexes = np.in1d(ids, platform_markers)

            if (np.sum(platform_idexes) == 4):
                rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
                B = np.average(tvecs[platform_idexes == False], axis=0)
                T = np.average(tvecs[platform_idexes == True], axis=0)
                if (ts.size == 0):
                    ts = B
                else:
                    ts = np.vstack((ts, B))

    print(ts.shape)
    print(np.std(ts[:,0]))
    print(np.std(ts[:,1]))
    print(np.std(ts[:,2]))

def main():
    extract_points_from_video()

if __name__ == '__main__':
    main()
