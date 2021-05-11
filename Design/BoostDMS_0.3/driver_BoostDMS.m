% driver_BoostDMS.m script file
%
% Purpose:
%
% File driver_BoostDMS applies the boostDSM algorithm to compute the 
% complete Pareto front for the bound constrained optimization problem
% ZDT1,described in E. Zitzler, K. Deb, and L. Thiele, "Comparison of 
% Multiobjective Evolutionary Algorithms: Empirical Results", 
% Evolutionary Computation 8(2): 173-195, 2000.
%
% The optimizer uses the default options specified in the file 
% parameters_BoostDMS.m. An output report is produced, both at the screen 
% and in the text file BoostDMS_report.txt (stored at the BoostDMS 
% directory).
%
% BoostDMS Version 0.1
% Copyright (C) 2019 C. P. Br???s and A. L. Cust???dio.
% 
% BoostDMS Version 0.2
% Copyright (C) 2020 C. P. Br???s, A. L. Cust???dio, V. M. Duarte, P. D. Medeiros, S. Tavares
% http://ferrari.dmat.fct.unl.pt/personal/alcustodio/BoostDFO.htm
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA,
% or see <https://www.gnu.org/licenses/>.
%
format compact;
% Name of the problem to generate the output files.
problemName = "problem1";
% Name of the folder to store the results. Should be in the main directory
foldername = "results";
% Create function handles. Empty arrays ( '[]' ) can be used instead of 
% handles, e.g. for unconstrained optimization. All files should be on the
% problems' folder, as well as the initial file ('init_file' in this
% example)
%cd("problems")
func_F = str2func('func_F');
func_C = str2func('func_C');
grad_C = str2func('grad_C');
%cd("..");
%
% To use the parallel version, set the option 'parallel_version' to 1 
% in the file 'parameters_BoostDMS.m' and uncomment the lines starting with
% '%%%'. 'numProcessors' should be replaced with the number of processors
% to use. Alternatively, 'parpool()' uses all available processors.
%
% Open the parallel pool
parpool(6);
%
[info,time,iter,iter_suc,list_size,func_eval] = BoostDMS(func_F,'init_file',func_C,grad_C,problemName,foldername);
%
% Close the parallel pool
delete(gcp('nocreate'));
%
% End of driver_BoostDMS.
%