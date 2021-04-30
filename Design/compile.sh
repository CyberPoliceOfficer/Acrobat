export OMP_NUM_THREADS=6
g++ -fopenmp evaluate_workspace.cpp -o evaluate_workspace -Wall
g++ -shared -o libworkspace.so -fPIC evaluate_workspace.cpp -fopenmp