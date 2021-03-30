import numpy as np
from scipy.spatial.transform import Rotation as R
import numpy.matlib

class RSS_cc:
    def __init__ (self, r_b = 0.5, r_m = 0.3, d_b = 0.2, d_m = 0.4, d = 0.7, h = 0.3, phi = 0.3491):
        # Physical parameters of the platform in meters and radians.
        
        self.r_b = r_b #Radious of the base
        self.r_m = r_m #Radious of the platform
        
        self.d_b = d_b #Lenght between base's anchors
        self.d_m = d_m #Lenght between platform's anchors
        
        self.d = d #Platform's arm lenght
        self.h = h #Servo's arm lenght
        
        self.phi = phi #Angle between servo's arm and platform's base
         
        # Compute vector bk and mk
        
        k = np.arange(1,7)
        n = np.floor((k-1)/2)
        
        theta_b = n*(2/3)*np.pi + np.power(-1,k)*np.arcsin(d_b/(2*r_b))
        theta_m = n*(2/3)*np.pi + np.power(-1,k)*np.arcsin(d_m/(2*r_m))
        
        self.b_k = [r_b * np.cos(theta_b), r_b * np.sin(theta_b), np.zeros(6)]
        self.m_k = [r_m * np.cos(theta_m), r_m * np.sin(theta_m), np.zeros(6)]
            
        # Compute beta and gamma

        self.beta_k = n*(2/3)*np.pi + np.power(-1,k)*np.pi/2
        self.phi_k = np.power(-1,k+1)*phi
        
    def inverse_kinematics (self, pose):
        Rot = R.from_euler('zyz', pose[3:], degrees=False)
        
        i_k = np.matlib.repmat(pose[:3], 6, 1).T + np.matmul(Rot.as_matrix(), self.m_k) - self.b_k
        
        f_k = (np.cos(self.beta_k)*i_k[0,:] + np.sin(self.beta_k)*i_k[1,:])*2*self.h

        e_k = (np.sin(self.beta_k)*np.sin(self.phi_k)*i_k[0,:] - np.cos(self.beta_k)*np.sin(self.phi_k)*i_k[1,:] + np.cos(self.phi_k)*i_k[2,:])*2*self.h
    
        g_k = np.power(np.linalg.norm(i_k, axis=0), 2) -(self.d**2 - self.h**2)
    
        return np.arcsin(g_k/np.sqrt(e_k**2 + f_k**2)) - np.arctan2(f_k,e_k)
    
    def inverse_kinematic_jacobian (self, pose, alpha_k):
        Rot = R.from_euler('zyz', pose[3:], degrees=False)
        T = np.matlib.repmat(pose[:3], 6, 1).T
        M_k = np.matmul(Rot.as_matrix(),self.m_k) + T
        h_k = self.h * np.array([np.sin(self.beta_k)*np.sin(self.phi_k)*np.sin(alpha_k) + np.cos(self.beta_k)*np.cos(alpha_k),
                        -np.cos(self.beta_k)*np.sin(self.phi_k)*np.sin(alpha_k) + np.sin(self.beta_k)*np.cos(alpha_k),
                        np.cos(self.phi_k)*np.sin(alpha_k)])
        H_k = self.b_k + h_k
        Jx = np.concatenate(((M_k-H_k).T, np.cross(M_k-T, M_k-H_k, axis=0).T), axis=1)
        u_k = -np.array([np.cos(self.beta_k+np.pi/2)*np.cos(self.phi_k),
              np.sin(self.beta_k+np.pi/2)*np.cos(self.phi_k),
              np.sin(self.phi_k)])

        # dot product
        dot = np.identity(6)/np.sum(np.cross(H_k-self.b_k, M_k - H_k, axis=0)*u_k, axis=0)
        return np.matmul(Jx, dot)

            