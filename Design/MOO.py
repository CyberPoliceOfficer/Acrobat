import ctypes
import pathlib
import math 
import numpy as np
import numpy.ctypeslib as ctl
import time
import mobopt as mo

class ObjectiveFunction:
    def __init__ (self, N_x = 50, N_y = 50, N_z = 50, joint_angle = 1, dim = 3, technique = 'discretization', c_integration = 'ctypes'):
        self.config = np.array([N_x, N_y, N_z, dim, np.pi/3], dtype=np.float32)
        if c_integration == 'CFFI':
            teste = 1
        else:
            # Load the shared library into ctypes
            libname = pathlib.Path().absolute() / "libworkspace.so"
            libworkspace = ctypes.CDLL(libname)

            self.evaluate_workspace = libworkspace.evaluate_workspace
            self.evaluate_workspace.argtypes = [ctypes.c_int, ctl.ndpointer(np.float32, flags='aligned, c_contiguous'), 
                                                        ctl.ndpointer(np.float32, flags='aligned, c_contiguous'), 
                                                        ctl.ndpointer(np.float32, flags='aligned, c_contiguous')]
        
        if technique == 'tree':
            # TBD, not implemented yet
            self.technique = 1
        else:
            self.technique = 0

    def evaluate (self, x):
        result = np.array([0, 0, 0, 0, 0], dtype=np.float32) #V, GTCI, GRCI, GTSI, GRSI
        self.evaluate_workspace(0, self.config, x, result)
        return result

if __name__ == "__main__":
    library = "MOBOpt"
    approach = "Local"
    NObj = 5
    Nparam = 8
    NIter = 3600
    bounds = np.zeros((Nparam,2), dtype=np.float32)

    if approach == "Local":
        # r_b
        bounds[0,0] = 0
        bounds[0,1] = 0
        # r_m
        bounds[1,0] = 0
        bounds[1,1] = 0
        # d_b
        bounds[2,0] = 0
        bounds[2,1] = 0
        # d_m
        bounds[3,0] = 0
        bounds[3,1] = 0
        # d
        bounds[4,0] = 0
        bounds[4,1] = 0
        # h
        bounds[5,0] = 0
        bounds[5,1] = 0
        # phi
        bounds[6,0] = 0
        bounds[6,1] = 0
        # beta
        bounds[7,0] = 0
        bounds[7,1] = 0


    # Test the objective function
    print("Objective function performance test")
    params = np.array([0.5, 0.3, 0.2, 0.6667, 0.7, 0.3, 0.3491, np.pi/2], dtype=np.float32) #r_b, r_m, d_b, d_m, d, h, phi, beta
    ObjFun = ObjectiveFunction()   
    t = time.time()
    result = ObjFun.evaluate(params)
    print('Time = ', time.time() - t)
    print('Result = ', result)

    print('Using ', library, ' for optimization')
    if library == "MOBOpt":
        Optimizer = mo.MOBayesianOpt(target=ObjFun.evaluate,
                                    NObj=NObj,
                                    pbounds=bounds)
        Optimizer.initialize(init_points=60)                            
        front, pop = Optimizer.maximize(n_iter=NIter)

