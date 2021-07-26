import numpy as np
import cv2 as cv
import glob
from scipy.spatial.transform import Rotation as R
import scipy
from scipy import optimize
from scipy.optimize import Bounds
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from RSS_cc import RSS_cc
from random import uniform

def extract_points_from_video ():
    draw = False
    aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
    aruco_side = 0.026
    n_ids = 11
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
    x_axis = np.array([26,27,28,29], dtype=np.int32)
    y_axis = np.array([26,21,16,11], dtype=np.int32)

    Rot = []
    Obv = []
    U = []

    for fname in videos:
        print(fname)
        npzfile = np.load(fname)
        images = npzfile['images']
        us = npzfile['us']
        for img, us in zip(images, us):

            # Detect aruco markers
            #img = cv.undistort(img, mtx, dist)
            corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img, aruco_dict, mtx, dist)
        
            
            platform_idexes = np.in1d(ids, platform_markers)

            if (ids.size == n_ids):
                rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
                tvecs = np.squeeze(tvecs)
                
                # Find axis
                x_base = np.squeeze(tvecs[np.in1d(ids, x_axis[-1])] - tvecs[np.in1d(ids, x_axis[0])])
                x_base = x_base/np.linalg.norm(x_base)
                y_base = np.squeeze(tvecs[np.in1d(ids, y_axis[-1])] - tvecs[np.in1d(ids, y_axis[0])])
                y_base = y_base/np.linalg.norm(y_base)
                z_base = np.cross(x_base, y_base)
                #print("Test base (should be ~0):", np.matmul(y_base,x_base.T))
                
                O = np.array([x_base,y_base,z_base]).T
                #print(O.T-np.linalg.inv(O))              
                
                x_platform = O @ np.squeeze(tvecs[np.in1d(ids, 13)] - tvecs[np.in1d(ids, 19)])
                x_platform = x_platform/np.linalg.norm(x_platform)
                y_platform = O @ np.squeeze(tvecs[np.in1d(ids, 14)] - tvecs[np.in1d(ids, 19)])
                y_platform = y_platform/np.linalg.norm(y_platform)
                z_platform = np.cross(x_platform, y_platform)
                #print("Test plat (should be ~0):", np.matmul(y_platform,x_platform.T))

                M = np.array([x_platform,y_platform,z_platform]).T

                # might be wrong, need to doublecheck later
                Rot.append(M)

                obv_tmp = O @ np.average(tvecs[platform_idexes == True], axis=0)
                ref_tmp = O @ np.squeeze(tvecs[np.in1d(ids, x_axis[0])].T)

                Obv.append(obv_tmp - ref_tmp)
                U.append(us)
                

            if draw:
                

                p = obv_tmp - ref_tmp

                fig = plt.figure()
                ax = Axes3D(fig)
                ax.set_xlim3d(-0.4, 0.4)
                ax.set_ylim3d(-0.4, 0.4)
                ax.set_zlim3d(-0.4, 0.4)


                ax.plot3D([0, 0.1], [0, 0], [0, 0], 'r')
                ax.plot3D([0, 0], [0, 0.1], [0, 0], 'g')
                ax.plot3D([0, 0], [0, 0], [0, 0.1], 'b')

                '''
                ax.plot3D([p[0], p[0] + 0.1*x_platform[0]], [p[1], p[1] + 0.1*x_platform[1]], [p[2], p[2] + 0.1*x_platform[2]], 'r')
                ax.plot3D([p[0], p[0] + 0.1*y_platform[0]], [p[1], p[1] + 0.1*y_platform[1]], [p[2], p[2] + 0.1*y_platform[2]], 'g')
                ax.plot3D([p[0], p[0] + 0.1*z_platform[0]], [p[1], p[1] + 0.1*z_platform[1]], [p[2], p[2] + 0.1*z_platform[2]], 'b')
                '''
                ux = M @ np.array([1,0,0])
                uy = M @ np.array([0,1,0])
                uz = M @ np.array([0,0,1])

                ax.plot3D([p[0], p[0] + 0.1*ux[0]], [p[1], p[1] + 0.1*ux[1]], [p[2], p[2] + 0.1*ux[2]], 'r')
                ax.plot3D([p[0], p[0] + 0.1*uy[0]], [p[1], p[1] + 0.1*uy[1]], [p[2], p[2] + 0.1*uy[2]], 'g')
                ax.plot3D([p[0], p[0] + 0.1*uz[0]], [p[1], p[1] + 0.1*uz[1]], [p[2], p[2] + 0.1*uz[2]], 'b')
                

                

                plt.show()  

                print(obv_tmp - ref_tmp)
                print(us)

                # Detect aruco markers
                rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)

                if (len(ids) > 0):
                    cv.aruco.drawDetectedMarkers(img, corners, ids)
                    rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
                    for i in range(len(ids)):
                        img = cv.aruco.drawAxis(img, mtx, dist, rvecs[i], tvecs[i],  0.05)
                cv.imshow('img', img)
                cv.waitKey()
    return np.asarray(Obv), np.asarray(Rot), np.asarray(U)

def generate_test_points(x_prototype, T_base):
    #x_prototype = np.array([0.05257, 0.04814, 0.1624, 0.0831, 0.3491, np.pi/2, 0.1175, 0.027])
    Manipulator = RSS_cc(r_b = x_prototype[0], r_m = x_prototype[1], d_b = x_prototype[2], d_m = x_prototype[3], d = x_prototype[6], h = x_prototype[7])
    u0 = np.array([1555, 1450, 1550, 1520, 1450, 1420])
    gamma = 2000/np.pi
    k = np.arange(1,7)
    num_points = 2000
    Rot = []
    Obv = []
    U = []
    corner = 0.04
    for iterator in range(num_points):
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
            if not isNaN:
                u = u0 + np.power(-1,k) * gamma*thetas
                rot = R.from_euler('zyz', p[3:], degrees=False).as_matrix()
                U.append(u)
                Rot.append(rot)
                Obv.append(p[0:3] + T_base)

    return np.asarray(Obv), np.asarray(Rot), np.asarray(U)

class Objective_Function:
    def __init__ (self, Obv, Rot, U):
        self.Obv = Obv
        self.R = Rot
        self.U = U
        self.n_obs = np.shape(U)[0]
        self.proto = None

    def evaluate(self, x):
        if np.size(x) == 21:
            return self.evaluate_1(x)
        if np.size(x) == (13):
            new_x = np.zeros(21)
            new_x[0:6] = x[0:6]
            '''
            new_x[6] = self.r_b
            new_x[7] = self.r_m
            new_x[8] = self.d_b
            new_x[9] = self.d_m
            new_x[10] = self.phi_o
            new_x[11] = self.beta_0
            new_x[12] = self.d
            new_x[13] = self.h
            '''
            new_x[6:14] = self.proto
            new_x[14:17] = x[6:9]
            new_x[17] = x[9]
            new_x[18:21] = x[10:13]
            return self.evaluate_1(new_x)



    def evaluate_1(self, x):
        # Input array u_0, u_1, u_2, u_3, u_4, u_5, r_b, r_m, d_b, d_m, phi_0, beta_0, d, h, T, theta_z
        u_0 = x[0:6]
        r_b = x[6]
        r_m = x[7]
        d_b = x[8]
        d_m = x[9]
        phi_0 = x[10]
        beta_0 = x[11]
        d = x[12]
        h = x[13]
        T_base = x[14:17]
        theta_z = x[17]
        T_platform = x[18:21]

        # Manipulator params
        gamma = np.pi/2000 #(Delta alpha)/(Delta u)
        k = np.arange(1,7)
        n = np.floor((k-1)/2)

        alpha_k = np.power(-1,k)*gamma*(self.U - u_0)

        beta_k = n*(2/3)*np.pi + np.power(-1,k)*beta_0
        phi_k = np.power(-1,k+1)*phi_0

        theta_b = n*(2/3)*np.pi + np.power(-1,k)*d_b
        theta_m = n*(2/3)*np.pi + np.power(-1,k)*d_m

        b_k = np.array([r_b * np.cos(theta_b), r_b * np.sin(theta_b), np.zeros(6)])
        m_k = np.array([r_m * np.cos(theta_m), r_m * np.sin(theta_m), np.zeros(6)])

        h_k = h * np.array([np.sin(beta_k)*np.sin(phi_k)*np.sin(alpha_k) + np.cos(beta_k)*np.cos(alpha_k),
                            -np.cos(beta_k)*np.sin(phi_k)*np.sin(alpha_k) + np.sin(beta_k)*np.cos(alpha_k),
                            np.cos(phi_k)*np.sin(alpha_k)])
        h_k = np.moveaxis(h_k, 0, 1)
        H_k = b_k.reshape((1,b_k.shape[0],b_k.shape[1])) + h_k


        #test = np.einsum('mij,jk->mik', self.R, R.from_euler('z', theta_z, degrees=False).as_matrix())
        Rm = R.from_euler('z', theta_z, degrees=False).as_matrix() @ self.R

        d_n = (self.Obv - T_base - Rm @ T_platform).reshape((self.Obv.shape[0],self.Obv.shape[1],1)) + Rm @ m_k - H_k
        
        return np.sum(((np.linalg.norm(d_n, axis=1) - d)**2))
    
        
def main():
    load_data = 2 # 0 - Load from memory; 1 - generate from images; 2 - generate from kinematics for testing
    opt_method = 2 # 1 - reduced number of variables ; 2 - All variables
    x_prototype = np.array([0.05257, 0.04814, 0.1624, 0.0831, 0.3491, np.pi/2, 0.1175, 0.027])
    if load_data == 1:
        Obv, Rot, U = extract_points_from_video()
        np.savez('TRB_data', Obv=Obv, Rot=Rot, U=U)
    
    if load_data == 0:
        npzfile = np.load('TRB_data.npz')
        Obv = npzfile['Obv']
        Rot = npzfile['Rot']
        U = npzfile['U']
    
    if load_data == 2:
        T_base =  np.array([0.25, 0.0, 0.04])
        Obv, Rot, U = generate_test_points(x_prototype, T_base)


    #print("Obv shape:",np.shape(Obv))
    #print("Rot shape:", np.shape(Rot))
    #print("U shape:",np.shape(U))
    #print(Obv)
    obj = Objective_Function(Obv, Rot, U)

    #(r_b = 0.05257, r_m = 0.04814, d_b = 0.1624, d_m = 0.0831, d = 0.1175, h = 0.027, phi = 0.3491, beta = np.pi/2
    # u_0, u_1, u_2, u_3, u_4, u_5, 
    # r_b, r_m, d_b, d_m, phi_0, beta_0, 
    # d, h, T_base, theta_z9
    # T_platform

    if opt_method == 1:
        x0 = np.array([1500, 1500, 1500, 1500, 1500, 1500, 
                        0.25, 0, 0, 0, 0, 0, 0])
        obj.proto = x_prototype
        bounds = Bounds(np.asarray([1200.0, 1200.0, 1200.0, 1200.0, 1200.0, 1200.0,
                         0.0, -0.05, 0, -np.pi/4,
                        -0.05, -0.05, -0.05]), 
                        np.asarray([1800.0, 1800.0, 1800.0, 1800.0, 1800.0, 1800.0, 
                        0.3, 0.05, 0.1, np.pi/4,
                        0.05, 0.05, 0.05]))
    if opt_method == 2:
        x0 = np.array([1500, 1500, 1500, 1500, 1500, 1500, 
        0.05257, 0.04814, 0.1624, 0.0831, 0.3491, np.pi/2,
        0.1175, 0.027, 0.25, 0, 0.03, 0, 
        0, 0, 0])
        bounds = Bounds(np.asarray([1300.0, 1300.0, 1300.0, 1300.0, 1300.0, 1300.0,
                        0.04, 0.04, 0, 0, 0.2, 1.55, 
                        0.09, 0.01, 0.05, -0.05, 0, -np.pi/4,
                        -0.05, -0.05, -0.05]), 
                        np.asarray([1700.0, 1700.0, 1700.0, 1700.0, 1700.0, 1700.0, 
                        0.06, 0.06, np.pi/3, np.pi/3, 0.6, 1.6, 
                        0.15, 0.04, 0.3, 0.05, 0.1, np.pi/4,
                        0.05, 0.05, 0.05]))

    print("F(x0) = ", obj.evaluate(x0))
    
    #OptimizeResult = scipy.optimize.minimize(obj.evaluate_1, x0, args=(), method='Nelder-Mead', tol=None, callback=None, options={'maxiter': 1e6, 'maxfev': None, 'disp': True, 'adaptive': True})
    OptimizeResult = scipy.optimize.minimize(obj.evaluate, x0, args=(), method='Powell', bounds=bounds, tol=None, callback=None, options={'xtol': 0.0001, 'ftol': 0.0001, 'maxiter': 1e6, 'maxfev': None, 'disp': True, 'direc': None, 'return_all': False})
    xp = OptimizeResult.x
    np.savez('OptimizeResult', OptimizeResult=OptimizeResult)
    if opt_method == 1:
        x = np.zeros(21)
        x[0:6] = xp[0:6]
        x[6:14] = x_prototype
        x[14:17] = xp[6:9]
        x[17] = xp[9]
        x[18:21] = xp[10:13]    
    else:
        x = xp

    print('u_0 = ', x[0:6])
    print('r_b = ', x[6])
    print('r_m = ', x[7])
    print('d_b = ', x[8])
    print('d_m = ', x[9])
    print('phi_0 = ', x[10])
    print('beta_0 = ', x[11])
    print('d = ', x[12])
    print('h = ', x[13])
    print('T_base = ', x[14:17])
    print('theta_z = ', x[17])
    print('T_mobile = ', x[18:21])
    #print((obj.Obv - obj.Ref -  x[14:17] ).reshape((obj.Obv.shape[0],obj.Obv.shape[1],1)))
if __name__ == '__main__':
    main()

