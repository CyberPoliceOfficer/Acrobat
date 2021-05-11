% parameters_BoostDMS.m script file
%
% Purpose:
%
% File parameters_BoostDMS sets the algorithmic strategies, parameters,
% constants, and tolerances values to be used by the function BoostDMS, 
% which can be user modified.
%
% DMS Version 0.2.
% Copyright (C) 2011 A. L. Cust???dio, J. F. A. Madeira, A. I. F. Vaz, 
% and L. N. Vicente.
% http://www.mat.uc.pt/dms
%
% BoostDMS Version 0.3.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stopping Criteria.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
stop_alfa = 1;      % 0-1 variable: 1 if the stopping criterion is based
                    % on the step size parameter; 0 otherwise.
tol_stop  = 10^-3;  % Lowest value allowed for the step size parameter.
%
stop_feval = 1;     % 0-1 variable: 1 if the stopping criterion is based
                    % on a maximum number of function evaluations; 0
                    % otherwise.
max_fevals = 40000;  % Maximum number of function evaluations allowed.
%
stop_fparcycles = 1;     % 0-1 variable: 1 if the stopping criterion is based
                         % on a maximum number of parallel function evaluation
                         % cycles; 0 otherwise. Can only be used with the 
                         % parallel version
max_fparCycles = 40000;  % Maximum number of parallel function 
                         %   evaluation cycles allowed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithmic Options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialization.
%
list = 3;   % 0-4 variable: 0 if the algorithm initializes the iterate
            % list with a single point; 1 if a latin hypercube sampling
            % strategy is considered for initialization; 2 if 
            % random sampling is used; 3 if points are considered 
            % equally spaced in the line segment, joining the 
            % variable upper and lower bounds and adding the centroid of
            % the feasible region; 4 if the algorithm is initialized with
            % a list provided by the optimizer;
%
user_list_size = 0;  % 0-1 variable: 1 if the user sets below the iterate 
                     % list initial size; 0 if the iterate list initial size
                     % equals the problem dimension, plus the centroid.
%                     
nPini = 30; % Number of points to consider in the iterate
            % list initialization, when its size is defined 
            % by the user.            
%
% Cache Use.
%
cache = 1; % 0-1 variable: 0 if point evaluation is always done;
           % 1 if a cache is maintained.
%
tol_match = tol_stop; % Tolerance used in point comparisons, when 
                      % considering a cache.
%
% Search step.
%
search_option = 1;  % 0-1 variable: 0 for no search step; 1 for a search step 
                    % based on quadratic polynomial models (cache should be
                    % set equal to 1).
%
regopt        = 1;  % 0-1 variable: 1 if, at the search step, a regression 
                    % model is considered when the number of points in the
                    % interpolation set exceeds the number of points 
                    % required for complete quadratic interpolation; 0 if 
                    % some of the interpolation points are discarded, in 
                    % order to only compute determined quadratic 
                    % interpolation models.
                    
all_subproblems = 0; % 0-1 variable: 0 computes the models level by level.
                     % 1 builds all quadratic models at the same time,
                     % after the first level. Only relevant if the problem
                     % has at least 3 objective functions
%                    
% Centre selection options.
%
spread_option = 1;      % 0-1 variable: 0 if for each point in the current approximation
                        % to the Pareto front, each component of the objective 
                        % function is projected in the corresponding dimension and 
                        % points are ordered according to the largest gap between 
                        % consecutive points, just before polling (standard gamma metric);
						% 1 if the largest gaps between consecutive points are
                        % ordered by objective function component, rotating
                        % on every iteration (cyclic gamma metric).
%
selection_strategy = 1; % 0-2 variable: 0 to select only one centre point
						% in each iteration, chosen by decreasing
                        % order of spread between points. 1 to select the
                        % first two points based on spread (biggest gap)
                        % and then one point with minimum value for each
                        % objective function component (total 2 + nobj centres). 
						% 2 same as 1, but new points obtained in the search
                        % step will be selected as poll centres, in
                        % addition to the centres already selected for
                        % polling.
%
% Directions and step size parameter.
%
dir_dense = 0;    % 0-1 variable: 1 if a dense set of directions should be 
                  % considered for polling; 0 otherwise.
%
alfa_ini  = 1;    % Initial step size. 
%
beta_par  = 0.5;  % Coefficient for step size contraction.
%
gamma_par = 1;    % Coefficient for step size expansion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel Options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
parallel_version = 1;   % 0-1 variable: 1 for the parallel version, 
                        % 0 for the serial one. This parallel version
                        % includes parallel objective function evaluations.
%
% End of parameters_BoostDMS.