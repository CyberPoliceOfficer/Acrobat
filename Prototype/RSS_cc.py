import numpy as np
from scipy.spatial.transform import Rotation as R
import numpy.matlib
import time

def divide_box (box, combination_key, dim = 3):
    boxes = np.zeros((2**dim,dim*2))
    i = 0
    for subbox in boxes:
        indices = format(i, combination_key)
        for j in range(dim):
            if (int(indices[j])):
                subbox[2*j] = box[2*j] + (box[2*j+1] - box[2*j])/2
                subbox[2*j+1] = box[2*j+1] 
            else:
                subbox[2*j] = box[2*j] 
                subbox[2*j+1] = box[2*j+1] - (box[2*j+1] - box[2*j])/2
        i = i + 1
    return boxes

def estimate_box_volume (box, dim = 3):
    volume = 1
    for i in range(dim):
        volume = volume*(box[i*2+1]-box[i*2])
    return volume

class RSS_cc:
    def __init__ (self, r_b = 0.5, r_m = 0.3, d_b = 0.2, d_m = 0.4, d = 0.7, h = 0.3, phi = 0.3491, beta = np.pi/2):
        # Physical parameters of the platform in meters and radians.       
        self.r_b = r_b #Radious of the base
        self.r_m = r_m #Radious of the platform
        
        self.d_b = d_b #Lenght between base's anchors
        self.d_m = d_m #Lenght between platform's anchors
        
        self.d = d #Platform's arm lenght
        self.h = h #Servo's arm lenght
        
        self.phi = phi #Angle between servo's arm and platform's base
        self.beta = beta
         
        # Compute vector bk and mk 
        k = np.arange(1,7)
        n = np.floor((k-1)/2)
        
        self.theta_b = n*(2/3)*np.pi + np.power(-1,k)*np.arcsin(d_b/(2*r_b))
        self.theta_m = n*(2/3)*np.pi + np.power(-1,k)*np.arcsin(d_m/(2*r_m))
        
        self.b_k = [r_b * np.cos(self.theta_b), r_b * np.sin(self.theta_b), np.zeros(6)]
        self.m_k = [r_m * np.cos(self.theta_m), r_m * np.sin(self.theta_m), np.zeros(6)]
            
        # Compute beta and gamma
        self.beta_k = n*(2/3)*np.pi + np.power(-1,k)*self.beta
        self.phi_k = np.power(-1,k+1)*phi
        
    def inverse_kinematics (self, pose):
        Rot = R.from_euler('zyz', pose[3:], degrees=False)
        
        i_k = np.matlib.repmat(pose[:3], 6, 1).T + np.matmul(Rot.as_matrix(), self.m_k) - self.b_k
        
        f_k = (np.cos(self.beta_k)*i_k[0,:] + np.sin(self.beta_k)*i_k[1,:])*2*self.h

        e_k = (np.sin(self.beta_k)*np.sin(self.phi_k)*i_k[0,:] - np.cos(self.beta_k)*np.sin(self.phi_k)*i_k[1,:] + np.cos(self.phi_k)*i_k[2,:])*2*self.h
    
        g_k = np.power(np.linalg.norm(i_k, axis=0), 2) -(self.d**2 - self.h**2)
    
        return np.arcsin(g_k/np.sqrt(e_k**2 + f_k**2)) - np.arctan2(f_k,e_k)
    
    def inverse_kinematic_jacobian (self, pose):
        Rot = R.from_euler('zyz', pose[3:], degrees=False)
        alpha_k = self.inverse_kinematics(pose)
        
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

    def home_pose (self):
        z = np.sqrt(self.d**2 - (self.r_m*np.cos(self.theta_m[0]) - self.r_b*np.cos(self.theta_b[0]) + self.h*np.cos(self.beta_k[0]))**2 - (self.r_m*np.sin(self.theta_m[0]) - self.r_b*np.sin(self.theta_b[0]) - self.h*np.sin(self.beta_k[0]))**2)         
        return np.array([0, 0, z, 0, 0, 0])
    
    def check_pose(self, pose):
        '''
            Checks if a given pose p is inside the workspace, considering kinematic and joint constraints.
            Returns True if Yes, False else
        '''
        # Kinematic constraints

        Rot = R.from_euler('zyz', pose[3:], degrees=False)
            
        i_k = np.matlib.repmat(pose[:3], 6, 1).T + np.matmul(Rot.as_matrix(), self.m_k) - self.b_k
        
        f_k = (np.cos(self.beta_k)*i_k[0,:] + np.sin(self.beta_k)*i_k[1,:])*2*self.h

        e_k = (np.sin(self.beta_k)*np.sin(self.phi_k)*i_k[0,:] - np.cos(self.beta_k)*np.sin(self.phi_k)*i_k[1,:] + np.cos(self.phi_k)*i_k[2,:])*2*self.h

        g_k = np.power(np.linalg.norm(i_k, axis=0), 2) - (self.d**2 - self.h**2)
    

        if(not np.all((np.abs(g_k) <= np.sqrt(e_k**2 + f_k**2)))):
            return False

    
        # Joint constraints.
        # TBD
        return True

    def check_bounded_box (self, box, combination_key, epsilon, dim = 3):
        '''
            Checks if a bounded box is inside the workspace. Returns True if Yes, False else.
            box must be either a 3D or 6D box.  
        '''
        
        
        
        vertexes = int(2**(dim)) # Either 8 or 64
        
        aux_pose = np.zeros(6)
        
        counter = vertexes
        
        for i in range(vertexes):
            indices = format(i, combination_key)
            for j in range(dim):
                aux_pose[j] = box [2*j + int(indices[j])]
            if (not self.check_pose(aux_pose)):
                counter = counter - 1
        
        # If the box is either completly out or completly inside
        if counter == 0:
            return 0
        if counter == vertexes:
            return estimate_box_volume (box, dim)
        
        # If the box volume is smaller than epsilon, don't search
        volume = estimate_box_volume (box, dim)
        if (volume < epsilon):
            return counter/vertexes*volume
        
        # If the box is in the workspace border, must be divided
        boxes = divide_box (box, combination_key, dim)
        
        volume = 0
        for subbox in boxes:
            subbox_volume = self.check_bounded_box (subbox, combination_key, epsilon, dim)
            volume = volume + subbox_volume
        return volume
            

    def estimate_workspace_metrics_tree (self, dim = 3, epsilon = 1e-5):
        '''
            Estimates the workspace and the workspace metrics inside it
        '''
        combination_key = "0" + str(dim) + "b"
        inner_point = self.home_pose()
        translation_maj = self.d + self.h - np.sin(np.pi/3)*self.r_b
        orientation_maj = np.pi/3
        
        # First boxes must be carefully assembled
        boxes = np.zeros((2**dim,dim*2))
        i = 0
        for subbox in boxes:
            indices = format(i, combination_key)
            for j in range(dim):
                if j < 2: # x and y
                    if (int(indices[j])):
                        subbox[2*j] = 0
                        subbox[2*j+1] = translation_maj
                    else:
                        subbox[2*j] = -1*translation_maj
                        subbox[2*j+1] = 0
                elif j == 2: # z
                    if (int(indices[j])):
                        subbox[2*j] = 0
                        subbox[2*j+1] = inner_point[2]
                    else:
                        subbox[2*j] = inner_point[2]
                        subbox[2*j+1] = self.d + self.h
                else: # 3 rotation angles
                    if (int(indices[j])):
                        subbox[2*j] = 0
                        subbox[2*j+1] = orientation_maj
                    else:
                        subbox[2*j] = -1*orientation_maj
                        subbox[2*j+1] = 0
            i = i + 1
        
        volume = 0
        for subbox in boxes:
            subbox_volume = self.check_bounded_box (subbox, combination_key, epsilon, dim)
            volume = volume + subbox_volume
        return volume

        
    def estimate_workspace_metrics_discretization (self, n_x = 100, n_y = 100, n_z = 50):
        '''
            Estimates the workspace and the workspace metrics inside it
        '''
        xy_maj = self.d + self.h - np.sin(np.pi/3)*self.r_b
        z_maj = self.d + self.h
        
        dx = 2*xy_maj/(n_x-1)
        dy = 2*xy_maj/(n_y-1)
        dz = z_maj/(n_z-1)
        
        print(dx*dy*dz)
        
        pose = np.zeros(6)
        volume = 0
        
        for x in range(n_x):
            for y in range(n_y):
                for z in range(n_z):
                    pose[0] = -1*xy_maj + x*dx
                    pose[1] = -1*xy_maj + y*dy
                    pose[2] = z*dz
                    if (self.check_pose(pose)):
                        volume = volume + dx*dy*dz
        
        return volume
                    

        
def main():
    t = time.time()
    manipulator = RSS_cc()
    volume = manipulator.estimate_workspace_metrics_tree(epsilon = 1e-07)
    #volume = manipulator.estimate_workspace_metrics_discretization()
    print('Volume:', volume)
    print('Time:', time.time() - t)



if __name__ == "__main__":
    main()