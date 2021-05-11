import ctypes
import pathlib
import math 
import numpy as np
import numpy.ctypeslib as ctl
import time
import mobopt as mo
from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_problem
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.factory import get_problem, get_sampling, get_crossover, get_mutation

class ObjectiveFunction:
    def __init__ (self, N_x = 50, N_y = 50, N_z = 50, joint_angle = 1, dim = 3, nFunObg = 3, technique = 'discretization', c_integration = 'ctypes', fixed_params = np.NaN*np.zeros(8)):
        self.config = np.array([N_x, N_y, N_z, dim, np.pi/3, nFunObg], dtype=np.float32)
        self.fixed_params = fixed_params
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
        fp = np.copy(self.fixed_params)
        fp[np.isnan(fp)]=x
        result = np.zeros(np.int(self.config[5]*2-1), dtype=np.float32) #V, GTCI, GRCI, GTSI, GRSI
        self.evaluate_workspace(0, self.config, fp.astype(np.float32), result)
        return result

if __name__ == "__main__":
    solver = "NSGA2"
    approach = "Normalized_d_h"
    NIter = 100

    # Test the objective function
    print("Objective function performance test")
    params = np.array([0.5, 0.3, 0.2, 0.6667, 0.7, 0.3, 0.3491, np.pi/2], dtype=np.float32) #r_b, r_m, d_b, d_m, d, h, phi, beta
    ObjFunTest = ObjectiveFunction()   
    t = time.time()
    result = ObjFunTest.evaluate(params)
    print('Time = ', time.time() - t)
    print('Result = ', result)

    print('Using ', solver, ' for optimization')

    # Define the solvers
    if solver == "MOBOpt":
        # Define the Problems
        if approach == "Normalized_d_h":
            NObj = 5
            Nparam = 2
            
            bounds = np.zeros((Nparam,2), dtype=np.float32)
            ObjFun = ObjectiveFunction(fixed_params=np.array([0.05257, 0.04814, 0.2, 0.6667, np.NaN, np.NaN, 0.3491, np.pi/2], dtype=np.float32)) 
            # d
            bounds[0,0] = 0.03
            bounds[0,1] = 0.2
            # h
            bounds[1,0] = 0.01
            bounds[1,1] = 0.03

        if approach == "Local":
            NObj = 5
            Nparam = 8
            bounds = np.zeros((Nparam,2), dtype=np.float32)
            ObjFun = ObjectiveFunction() 
            # r_b
            bounds[0,0] = 0.01
            bounds[0,1] = 0.06
            # r_m
            bounds[1,0] = 0.01
            bounds[1,1] = 0.06
            # d_b
            bounds[2,0] = 0
            bounds[2,1] = np.pi/3
            # d_m
            bounds[3,0] = 0
            bounds[3,1] = np.pi/3
            # d
            bounds[4,0] = 0.01
            bounds[4,1] = 0.2
            # h
            bounds[5,0] = 0.01
            bounds[5,1] = 0.2
            # phi
            bounds[6,0] = -np.pi/2
            bounds[6,1] = np.pi/2
            # beta
            bounds[7,0] = 0
            bounds[7,1] = np.pi
        
        Optimizer = mo.MOBayesianOpt(target=ObjFun.evaluate, NObj=NObj, pbounds=bounds, verbose=True)
        Optimizer.initialize(init_points=1000)
        print('Innit done, optimizing... ')                               
        front, pop = Optimizer.maximize(n_iter=NIter, FrontSampling=[50,100,200,400,800])
        np.savez(solver + '_' + approach, front=front, pop=pop)
    

    if solver == "NSGA2":
        n_gens = 500
        class OtimizeRSS(Problem):
            def __init__(self, n_var = 8, n_obj = 3, n_constr = 1, 
            xl=np.array([0.01, 0.01, 0, 0, 0.01, 0.01, -np.pi/2, 0]),
            xu=np.array([0.06, 0.06, np.pi/3, np.pi/3, 0.2, 0.2, np.pi/2, np.pi])):
                super().__init__(n_var=n_var,
                                n_obj=n_obj,
                                n_constr=n_constr,
                                xl=xl,
                                xu=xu)
                self.ObjFun = ObjectiveFunction(nFunObg = 2)
                self.n_obj = n_obj

            def _evaluate(self, X, out, *args, **kwargs):
                #r_b = X[:, 0], r_m = X[:, 1], d_b = X[:, 2], d_m = X[:, 3] , d = X[:, 4], h = X[:, 5], phi = X[:, 6], beta = X[:, 7]
                num_rows, num_cols = X.shape
                tmp = np.empty ((num_rows,self.n_obj), dtype=np.float32)
                for i in range(num_rows):
                    tmp[i,:] = self.ObjFun.evaluate (X[i,:])
                
                out["F"] = tmp
                out["G"] = -X[:, 4]**2 + (X[:, 1]*np.cos(-X[:, 3]) - X[:, 0]*np.cos(-X[:, 2]) + X[:, 5]*np.cos(X[:, 7]))**2 
                + (X[:, 1]*np.sin(-X[:, 3]) -X[:, 0]*np.sin(X[:, 2]) - X[:, 5]*np.sin(-X[:, 7]))**2
        

        # Get Algorithm
        problem = OtimizeRSS()
        algorithm = NSGA2(pop_size=500,
                  sampling=get_sampling("bin_random"),
                  crossover=get_crossover("bin_two_point"),
                  mutation=get_mutation("bin_bitflip"),
                  eliminate_duplicates=True)       

        np.save("checkpoint", algorithm)
        for i in range(0,n_gens,10):
            checkpoint, = np.load("checkpoint.npy", allow_pickle=True).flatten()
            print("Loaded Checkpoint:", checkpoint)
            res = minimize(problem, algorithm, ('n_gen', i),
                            seed=1, verbose=True, copy_algorithm=False)
            np.savez(solver + '_' + approach, res=res)
        

        plot = Scatter()
        plot.add(problem.pareto_front(), plot_type="line", color="black", alpha=0.7)
        plot.add(res.F, color="red")
        plot.show()
                        
