import numpy as np
import cv2 as cv
import glob
from scipy.spatial.transform import Rotation as R
import scipy
from scipy import optimize
from scipy.optimize import Bounds

def extract_points_from_video ():
    aruco_dict = cv.aruco.Dictionary_get(cv.aruco.DICT_6X6_1000)
    aruco_side = 0.026
    n_ids = 18
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
    x_axis = np.array([30,31,32,33], dtype=np.int32)
    y_axis = np.array([15,20,25,30], dtype=np.int32)

    R = []
    B = []
    T = []
    U = []

    for fname in videos:
        print(fname)
        npzfile = np.load(fname)
        images = npzfile['images']
        us = npzfile['us']
        for img, us in zip(images, us):

            # Detect aruco markers
            corners, ids, rejectedImgPoints	= cv.aruco.detectMarkers(img, aruco_dict, mtx, dist)
            
            platform_idexes = np.in1d(ids, platform_markers)

            if (ids.size == n_ids):
                rvecs, tvecs, _objPoints = cv.aruco.estimatePoseSingleMarkers(corners, aruco_side, mtx, dist)
                tvecs = np.squeeze(tvecs)
                
                # Find axis
                x_base = np.squeeze(tvecs[np.in1d(ids, 33)] - tvecs[np.in1d(ids, 30)])
                x_base = x_base/np.linalg.norm(x_base)
                y_base = np.squeeze(tvecs[np.in1d(ids, 15)] - tvecs[np.in1d(ids, 30)])
                y_base = y_base/np.linalg.norm(y_base)
                z_base = np.cross(y_base,x_base)
                #print("Test (should be ~0):", np.matmul(y_base,x_base.T))
                
                O = np.array([x_base,y_base,z_base])
                #print(O.T-np.linalg.inv(O))              
                
                x_platform = np.squeeze(tvecs[np.in1d(ids, 13)] - tvecs[np.in1d(ids, 19)])
                y_platform = np.squeeze(tvecs[np.in1d(ids, 14)] - tvecs[np.in1d(ids, 19)])
                z_platform = np.cross(y_platform,x_platform)
                # print("Test (should be ~0):", np.matmul(y_platform,x_platform.T))

                M = np.array([x_platform,y_platform,z_platform])

                # might be wrong, need to doublecheck later
                R.append(O.T @ M) 

                # Fit a vector for axis
                '''
                X_norm = tvecs[np.in1d(ids, x_axis)] - np.average(tvecs[np.in1d(ids, x_axis)])
                ux, s, v = np.linalg.svd((1./X_norm.shape[0])*np.matmul(X_norm.T,X_norm))
                x_u = ux[0,:]
                print(x_u.shape)

                Y_norm = tvecs[np.in1d(ids, y_axis)] - np.average(tvecs[np.in1d(ids, y_axis)])
                uy, s, v = np.linalg.svd((1./Y_norm.shape[0])*np.matmul(Y_norm.T,Y_norm))
                y_u = uy[0,:]
                print(y_u.shape)
                
                
                print(np.matmul(y_u,x_u))
                '''

                T.append(np.average(tvecs[platform_idexes == True], axis=0))
                B.append(np.average(tvecs[platform_idexes == False], axis=0))
                U.append(us)

    return T, R, B, U

class Objective_Function:
    def __init__ (self, T, R, B, U):
        self.T = T
        self.R = R
        self.B = B
        self.U = U
        self.n_obs = U.shape[0]

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

        # Manipulator params
        gamma = np.pi/2000 #(Delta alpha)/(Delta u)
        k = np.arange(1,7)
        n = np.floor((k-1)/2)

        alpha_k = gamma*(self.U - u_0)
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

        Rm = self.R @ R.from_euler('z', theta_z, degrees=False).as_matrix()

        d_n = (self.T + T_base).reshape((self.T.shape[0],self.T.shape[1],1)) + Rm @ m_k - H_k
        return np.sum(((np.linalg.norm(d_n, axis=1) - d)**2))
        


def main():
    load_data = False
    if load_data:
        T, R, B, U = extract_points_from_video()
        np.savez('TRB_data', T=T, R=R, B=B, U=U)
    
    npzfile = np.load('TRB_data.npz')
    T = npzfile['T']
    R = npzfile['R']
    B = npzfile['B']
    U = npzfile['U']
    
    obj = Objective_Function(T,R,B,U)

    obj.evaluate_1(np.ones(18))
    #(r_b = 0.05257, r_m = 0.04814, d_b = 0.1624, d_m = 0.0831, d = 0.1175, h = 0.027, phi = 0.3491, beta = np.pi/2
    # u_0, u_1, u_2, u_3, u_4, u_5, 
    # r_b, r_m, d_b, d_m, phi_0, beta_0, 
    # d, h, T, theta_z
    x0 = np.array([1500, 1500, 1500, 1500, 1500, 1500, 0.05257, 0.04814, 0.1624, 0.0831, 0.3491, np.pi/2,  0.1175, 0.027, 0, 0, 0, 0])
    bounds = Bounds([0, 0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 0, 0, 
                     0, 0, -2.0, -2.0, -2.0, -np.pi/4], 
                    [3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 
                    1.0, 1.0, np.pi/3, np.pi/3, np.pi/4, 3*np.pi/2, 
                    1.0, 1.0, -2.0, -2.0, -2.0, -np.pi/4])
    #OptimizeResult = scipy.optimize.minimize(obj.evaluate_1, x0, args=(), method='Nelder-Mead', tol=None, callback=None, options={'maxiter': 1e6, 'maxfev': None, 'disp': True, 'adaptive': True})
    OptimizeResult = scipy.optimize.minimize(obj.evaluate_1, x0, args=(), method='Powell', bounds=bounds, tol=None, callback=None, options={'xtol': 0.0001, 'ftol': 0.0001, 'maxiter': 1e6, 'maxfev': None, 'disp': True, 'direc': None, 'return_all': False})
    np.savez('OptimizeResult', OptimizeResult=OptimizeResult)
    print(OptimizeResult)

if __name__ == '__main__':
    main()
