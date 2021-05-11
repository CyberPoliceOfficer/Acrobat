    function [success_info,time,iter,iter_suc,list_size,func_eval,func_parCycles] = BoostDMS(func_F,file_ini,func_C,grad_C,problemName,foldername)
%
% Purpose:
%
% Function BoostDMS applies a direct multisearch method to the derivative-free
% multiobjective optimization problem:
%
%    min F(x) = (f_1(x),f_2(x),...,f_m(x))  s.t.  c_j(x) <= 0, j = 1,...,p,
%
% where x is a real vector of dimension n. The derivatives of the functions
% f_i, i = 1,..., m, are not used. Only function values are provided for F
% and C = (c_1,c_2,...,c_p).
%
% The user must provide: func_F (for F function values),
%
%                        and, if there are constrains other than bounds,
%
%                        func_C (for c_j, j = 1,...,p, function values)
%                        grad_C (for c_j, j = 1,...,p, gradient values).
%
% Input:
%
%         func_F (Function handle for the problem file.)
%
%         file_ini (Name of the file used for initialization.)
%
%		  func_C (Function handle for computing other problem constraints.)
%
%         grad_C (Function handle for computing gradients of problem constraints,
%                other than bounds.)
%
%		  problemName (Text to build BoostDMS_report.)
%
%		  foldername (Name of directory where files will be saved.)
%
% Output:
%
%         success_info (Success or error message.)
%
%         time (Computational time measured.)
%
%         iter (Number of iterations performed.)
%
%         iter_suc (Number of successful iterations performed.)
%
%         list_size (Number of points in the pareto front
%								at the end of the optimization.)
%
%         func_eval (Total number of function evaluations.)
%
%         func_parCycles (Total number of paralell function evaluation cycles.)
%
%         BoostDMS_report.txt, which records the iteration results.
%		  BoostDMS_plot.jpg, a plot of the retrieved Pareto front
%			for problems with 2 or 3 objective functions.
%		  BoostDMS_plot_matlab.fig, the full matlab figure corresponding
%           to the obtained plot.
%         BoostDMS_lastiteration.mat, which records all relevant
%				variables from the last iteration; this allows
%				to continue solving a problem at a later time
%
% Functions called: func_F, func_C, grad_C (application, user provided)
%
%                   match_point, const_model, obj_model, paretodominance,
%                   individual_models, sort_gamma (provided by the optimizer).
%
% DMS Version 0.3
% Copyright (C) 2015 A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz,
% and L. N. Vicente.
% http://www.mat.uc.pt/dms
%
% BoostDMS Version 0.2
% Copyright (C) 2020 C. P. Br�s, A. L. Cust�dio, V. M. Duarte, P. D. Medeiros, S. Tavares
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
%
% Try/catch last iteration variables (in case of unexpected error).
%
try
    time = clock;
    format long;
    warning off all;
    success_info = "success";
    iter = 0;
    list_size = 0;
    % Check if final time is already measured
    time_set = false;
    %
    % Load algorithmic strategies, parameters, constants, and tolerances
    % values.
    %
    parameters_BoostDMS;
    numCentres = 1;
    pollOnSucSearch = 0;
    if selection_strategy == 2
        pollOnSucSearch = 1;
    end
    %
    % Read the file used for initialization and initialize variables.
    % 	Bounds and at least one initial point will be included
    % 	if list==4 all information from the previous run
    % 		is uploaded.
    %
    if (list == 4)
        cd("problems");
        load(file_ini); %.mat file
        cd(fullfile("../"));
        %         x_initial = Plist;
        %         F_ini = Flist;
        %         alfa_list = alfa;
        max_fevals = max_fevals + func_eval;
        if parallel_version
            max_fparCycles = max_fparCycles + func_parCycles;
        end
        start_iter = iter;
        %
        % Algorithm counters and cache variables
        %	are implicitly loaded
        %
    else
        if cache
            CacheP     = [];
            CachenormP = [];
            CacheF     = [];
            Cache_infeasible = [];
            Cache_normInf = [];
        end
        start_iter = 0;
        iter_suc   = 0;
        %cd("problems");
        run(file_ini); %.m file
        %cd(fullfile("../"));
        lbound = lowerbound;
        ubound = upperbound;
        gamma_cycle = 1;
    end
    halt     = 0;
    % Find out the number of workers in the current parallel pool
    if parallel_version
        p = gcp('nocreate');
        numWorkers = p.NumWorkers;
    end
    %
    % Parameters used for model computation (not to be changed by user).
    %
    sigma      = 2;      % Coefficient used in trust-region definition.
    beta_trust = 3;      % Coefficient used in trust-region definition, only for
    % points selection.
    perc = 0.8;          % Percentage of nearest points to the current iterate,
    % which should be used in model computation.
    tol_Delta  = 10^-5;  % Minimum trust-region radius used in model computation.
    if search_option
        Opt = optimset('Display','off','GradConstr','on','GradObj','on','algorithm','sqp');
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization Step.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Define the problem size.
    %
    if ~isempty(lbound)
        n = size(lbound,1);
    else
        if ~isempty(ubound)
            n = size(ubound,1);
        else
            if ~isempty(x_initial)
                n = size(x_initial,1);
            else
                success_info = "An initial point or variable bounds should be provided.";
                return
            end
        end
    end
    %
    % Define the initial iterate list.
    %
    if (isempty(lbound) || isempty(ubound) || ((sum(isfinite(lbound))+...
            sum(isfinite(ubound)))~= (2*size(lbound,1)))) && isempty(x_initial)
        success_info = "An initial point or complete variable bounds should be provided.";
        return
    else
        if list == 0
            if ~isempty(x_initial)
                Pini = [x_initial];
            else
                Pini = (lbound + ubound)/2;
            end
        elseif list ~=4
            %
            %     Set seed for random strategies used in the iterate list initialization.
            %
            rand('state', sum(100*clock));
            %
            if user_list_size == 0
                nPini = n;
            end
            %
            if isempty(lbound)
                lbound_gen = x_initial - ones(n,1);
            else
                index      = isfinite(lbound);
                lbound_gen = lbound;
                if isfinite(x_initial)
                    lbound_gen(~index) = x_initial(~index) - max(x_initial(index)-lbound(index));
                end
            end
            if isempty(ubound)
                ubound_gen = x_initial + ones(n,1);
            else
                index      = isfinite(ubound);
                ubound_gen = ubound;
                if isfinite(x_initial)
                    ubound_gen(~index) = x_initial(~index) + max(ubound(index)-x_initial(index));
                end
            end
            if (list == 1)
                Pini = repmat(lbound_gen,1,nPini)+lhsdesign(nPini,n)'.*...
                    repmat((ubound_gen-lbound_gen),1,nPini);
            else
                if (list == 2)
                    Pini = repmat(lbound_gen,1,nPini)+rand(n,nPini).*...
                        repmat((ubound_gen-lbound_gen),1,nPini);
                else
                    if (nPini == 1)
                        Pini = (lbound_gen + ubound_gen)/2;
                    else
                        Pini = repmat(lbound_gen,1,nPini)+repmat([0:(nPini-1)]/(nPini-1),n,1)...
                            .*repmat((ubound_gen-lbound_gen),1,nPini);
                        if ~mod(nPini,2) && ~user_list_size
                            Pini = [Pini (lbound_gen + ubound_gen)/2];
                        end
                    end
                end
            end
        end
    end
    %
    if list ~= 4
        Flist     = [];
        Plist     = [];
        alfa      = [];
        func_eval = 0;
        func_parCycles = 0;
    end
    %
    % Evaluate the initial iterate list.
    %
    if list ~= 4
        mask_pte = true(1,size(Pini,2));
        if n > 1
            norm_pte = vecnorm(Pini,1);
        else
            norm_pte = abs(Pini);
        end
        for i=1:size(Pini,2)
            x_initial = Pini(:,i);
            %
            %  Check feasibility.
            %
            bound = [];
            if ~isempty(ubound)
                bound = x_initial - ubound;
            end
            if ~isempty(lbound)
                bound = [bound ; - x_initial + lbound];
            end
            if ~isempty(bound)
                m  = size(bound,1);
                if sum(bound <= 0) ~= m
                    mask_pte(i) = false;
                end
            end
            if mask_pte(i) && ~isempty(func_C)
                c_const = func_C(x_initial);
                m       = size(c_const,1);
                if sum(c_const <= 0) ~= m
                    mask_pte(i) = false;
                end
            end
            if mask_pte(i) && i < size(Pini,2)
                %
                % Verify if next points would match this
                % one in the cache, if they were to be
                % stored one by one (every point should be
                % at least 'tol_match' from eachother).
                %
                index = find(abs(norm_pte(i+1:end) - norm_pte(i)) <= tol_match) + i;
                if ~isempty(index)
                    auxP = Pini(:,index);
                    index2 = max(abs(auxP-x_initial),[],1) <= tol_match;
                    if any(index2)
                        index = index(index2);
                        mask_pte(index) = false;
                    end
                end
            end
        end
        pointsToEval = Pini(:,mask_pte);
        if ~all(mask_pte)
            norm_pte = norm_pte(mask_pte);
        end
        %
        % Evaluate all feasible points.
        %
        ftemp_pte(:,1) = func_F(pointsToEval(:,1));
        numEvalPoints = size(pointsToEval,2);
        if stop_feval && (func_eval + numEvalPoints > max_fevals)
            numEvalPoints = max_fevals - func_eval;
            halt = 1;
        end
        nobj = size(ftemp_pte,1);
        if numEvalPoints > 1
            ftemp_pte(:,numEvalPoints) = zeros(nobj, 1);
            if parallel_version && numEvalPoints > 2
                parfor (i = 2:numEvalPoints, numEvalPoints)
                    ftemp_pte(:,i) = func_F(pointsToEval(:,i)); %#ok<*PFBNS>
                end
            else
                for i = 2:numEvalPoints
                    ftemp_pte(:,i) = func_F(pointsToEval(:,i));
                end
            end
        end
        func_eval = func_eval + numEvalPoints;
        if parallel_version
            func_parCycles = func_parCycles + ceil((numEvalPoints-1)/numWorkers) + 1;
        end
        % Check if objective function values are finite
        mask_pte = sum(isfinite(ftemp_pte),1) == size(ftemp_pte,1);
        if halt || ~all(mask_pte)
            % Fill the infeasible points' cache
            if cache
                Cache_infeasible = pointsToEval(:,~mask_pte);
                Cache_normInf = norm_pte(~mask_pte);
            end
            % Crop results
            pointsToEval = pointsToEval(:,mask_pte);
            ftemp_pte = ftemp_pte(:,mask_pte);
            norm_pte = norm_pte(mask_pte);
        end
        if ~isempty(pointsToEval)
            if cache
                CacheP     = pointsToEval;
                CachenormP = norm_pte;
                CacheF     = ftemp_pte;
            end
            %
            % Check for nondominated points
            %
            Plist = pointsToEval(:,1);
            Flist = ftemp_pte(:,1);
            for i = 2:size(pointsToEval,2)
                [pdom,index_ndom] = paretodominance(ftemp_pte(:,i),Flist);
                if ~pdom
                    Plist = [Plist(:,index_ndom), pointsToEval(:,i)];
                    Flist = [Flist(:,index_ndom), ftemp_pte(:,i)];
                end
            end
            alfa = zeros(1,size(Plist,2)) + alfa_ini;
        end
    end
    %
    % Check if the iterate list is not empty.
    %
    if isempty(Flist)
        success_info = "The optimizer did not generate a feasible point or the initial point provided is not feasible.";
        return
    end
    %
    % Recalculate the number of centre points according to the selection
    % strategy
    %
    if selection_strategy == 1 || selection_strategy == 2
        numCentres = 2 + size(Flist,1);
    end
    %
    % Set seed for random generation of poll directions.
    %
    if dir_dense
        rand('state',1);
    end
    %
    % Clear used variables.
    %
    clear x_initial Pini pointsToEval ftemp_pte norm_pte mask_pte;
    %
    % Print the iteration report header.
    %
    fprintf('Iteration Report: \n\n');
    fprintf('| iter  | search_success | success | list size |     min alpha    |    max alpha     |\n');
    print_format = ['| %5d |       %2s       |    %2s   |   %5d   | %+13.8e  | %+13.8e  |\n'];
    fprintf(print_format, iter, '--','--', size(Plist,2), min(alfa), max(alfa));
    cd(foldername);
    output_file = "BoostDMS_report_" + problemName + ".txt";
    fresult = fopen(output_file,'w');
    cd(fullfile("../"));
    fprintf(fresult,'Iteration Report: \n\n');
    fprintf(fresult,'| iter  | search_success | success | list size |     min alpha    |     max alpha    |\n');
    fprintf(fresult,print_format, iter, '--','--',size(Plist,2), min(alfa), max(alfa));
    %
    if list ~= 4
        max_D  = -1; % max_D (Maximum size of the directions considered in the previous poll set.)
    end
    %
    while (~halt)
        func_iter = 0;
        move      = 0;
        success   = 0;
        search_success = 0;
        numNewCentres = 0;
        poll      = 1;
        %
        % Order the List of points based on the selected spread option
        %
        [Plist,Flist,alfa,gamma_cycle] = sort_gamma(Plist,Flist,alfa,stop_alfa,...
            tol_stop,spread_option,gamma_cycle);
        if (alfa(1) >= tol_stop) || ~stop_alfa
            orig_centres = Plist(:,1);
            centres = orig_centres;
            orig_alfa_centres = alfa(1);
            alfa_centres = orig_alfa_centres;
            numSelectedCentres = 1;
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     Search Step.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if max_D > 0 && search_option && cache
                %
                % Select all model centers according to the selection strategy
                %
                if selection_strategy == 1 || selection_strategy == 2
                    [orig_centres,orig_alfa_centres] = centre_selection(Plist,Flist,alfa,stop_alfa,tol_stop,numCentres);
                    numSelectedCentres = size(orig_centres,2);
                    cell_Paux = cell(numSelectedCentres,1);
                    cell_Faux = cell(numSelectedCentres,1);
                end
                Delta  = max(orig_alfa_centres * sigma * max_D, tol_Delta);
                orig_new_suc_alfa = orig_alfa_centres * gamma_par;
                for i = 1:numSelectedCentres
                    %
                    % Identify the current iterate in the cache.
                    %
                    index  = find(sum(CacheP == orig_centres(:,i),1) == size(Plist,1), 1);
                    P_aux  = [CacheP(:,index), CacheP(:,1:index-1), CacheP(:,index+1:size(CacheF,2))];
                    F_aux  = [CacheF(:,index), CacheF(:,1:index-1), CacheF(:,index+1:size(CacheF,2))];
                    %
                    % Order the cache by proximity to the current iterate (model
                    % center).
                    %
                    if size(CacheF,2) > (n+1)*(n+2)/2
                        [norm_aux,index] = sort(sum((P_aux-orig_centres(:,i)).^2,1),2,'ascend');
                        P_aux   = P_aux(:,index);
                        F_aux   = F_aux(:,index);
                        %
                        if regopt
                            p_max = (n+1)*(n+2);
                        else
                            p_max = (n+1)*(n+2)/2;
                        end
                        index  = norm_aux <= (beta_trust*Delta(i))^2;
                        P_aux1 = P_aux(:,index);
                        F_aux1 = F_aux(:,index);
                        P_aux1 = P_aux1(:,1:min(size(F_aux1,2),floor(perc*p_max)));
                        F_aux1 = F_aux1(:,1:min(size(F_aux1,2),floor(perc*p_max)));
                        %
                        index  = norm_aux > (beta_trust*Delta(i))^2;
                        P_aux2 = P_aux(:,index);
                        F_aux2 = F_aux(:,index);
                        P_aux  = [P_aux1, P_aux2(:,1:min(size(F_aux2,2),p_max-size(F_aux1,2)))];
                        F_aux  = [F_aux1, F_aux2(:,1:min(size(F_aux2,2),p_max-size(F_aux1,2)))];
                    end
                    if numSelectedCentres > 1
                        cell_Paux{i} = P_aux;
                        cell_Faux{i} = F_aux;
                    end
                end
                if numSelectedCentres > 1
                    P_aux = cell_Paux;
                    F_aux = cell_Faux;
                end
                %
                % Generate the quadratic models for each iterate and each
                % function component
                %
                [quad,P_search,cut_index,g_old,H_old] = individual_models(tol_match,Delta,P_aux,...
                    F_aux,orig_centres);
                %
                if quad
                    %
                    % Multiobjective minimization of the models, inside each level.
                    %
                    vmain = 1:nobj;
                    level = 1;            % Start at the minimum level
                    % Crop centres that didn't generate valid models
                    if ~all(cut_index)
                        centres = orig_centres(:,cut_index);
                        new_suc_alfa = orig_new_suc_alfa(cut_index);
                        Delta = Delta(cut_index);
                    else
                        centres = orig_centres;
                        new_suc_alfa = orig_new_suc_alfa;
                    end
                    numSelectedCentres = size(centres,2);
                    % Variable that tracks the centre that originated a
                    % given new point
                    centre_origin = repelem(1:numSelectedCentres,nobj);
                    % Variable to track success on each centre
                    successful_centres = false(1,numSelectedCentres);
                    if pollOnSucSearch
                        newPollCentres = [];
                        newPollCentres_alfa = [];
                    end
                    while ~halt && ~success && level <= nobj
                        %
                        %           Consider the quadratic models m_j(x+h), j=1,...,nobj defined in
                        %           individual models (where x is the center point) and solve the problem
                        %
                        %                   min  zeta
                        %                   s.t. m_j(x+h) <= zeta,  j in the set of combinations of nobj functions by level
                        %                       ||h||  <= Delta
                        %                       lbound <= x+h <= ubound
                        %
                        if level > 1
                            activeCentreIndexes = find(~successful_centres);
                            numSelectedCentres = size(activeCentreIndexes,2);
                            % Slice auxiliary variables to minimize data
                            % transfer to workers, before using fmincon
                            if level == 2
                                if iscell(F_aux)
                                    tmp_F_aux = zeros(nobj,numSelectedCentres);
                                    for j = 1:numSelectedCentres
                                        tmp_F_aux(:,j) = F_aux{activeCentreIndexes(j)}(:,1);
                                    end
                                    F_aux = tmp_F_aux;
                                else
                                    F_aux = F_aux(:,1);
                                end
                            end
                            if all_subproblems && nobj > 2
                                %
                                % Generate all combinations for all levels
                                % (above the first one).
                                %
                                vv = nchoosek(vmain,2);
                                cellvv = mat2cell(vv, ones(1, size(vv,1)));
                                for i = 3:nobj
                                    vv = nchoosek(vmain,i);
                                    cellvv = vertcat(cellvv, mat2cell(vv, ones(1, size(vv,1))));
                                end
                                level = nobj;
                            else
                                %
                                % All possible combinations of components of the objective
                                % function, for the level considered.
                                %
                                cellvv = nchoosek(vmain,level);
                            end
                            models_per_centre = size(cellvv,1);
                            numModelsEval = models_per_centre * numSelectedCentres;
                            P_search = zeros(size(Plist,1),numModelsEval);
                            mask_p = false(1,numModelsEval);
                            centre_origin = repelem(activeCentreIndexes,models_per_centre);
                            for k = 1:numModelsEval
                                centre_num = ceil(k/models_per_centre);
                                centre_idx = activeCentreIndexes(centre_num);
                                model_idx = mod(k-1,models_per_centre) + 1;
                                % Select one combination of the level considered.
                                if iscell(cellvv)
                                    v = cellvv{model_idx};
                                else
                                    v = cellvv(model_idx,:);
                                end
                                maxf  = max(F_aux(v,centre_num));   % Initial value for the variable zeta.
                                [p,~,exitflag] = fmincon(@(x)obj_model(x),[zeros(n,1);maxf],[],[],...
                                    [],[],[lbound-centres(:,centre_idx);-inf],[ubound-centres(:,centre_idx);maxf],...
                                    @(x)const_model(x,v,F_aux(:,centre_num),g_old(:,:,centre_idx),H_old(:,:,:,centre_idx),Delta(centre_idx),func_C,grad_C),Opt);
                                if exitflag == 1
                                    P_search(:,k) = centres(:,centre_idx) + p(1:n);
                                    mask_p(k) = true;
                                end
                            end
                            if ~all(mask_p)
                                P_search = P_search(:,mask_p);
                                centre_origin = centre_origin(mask_p);
                            end
                        end
                        if ~isempty(P_search)
                            P_search = min(max(lbound,P_search),ubound);
                            mask_pte = true(1,size(P_search,2));
                            norm_pte = vecnorm(P_search,1,1);
                            for i = 1:size(P_search,2)
                                if mask_pte(i)
                                    xtemp = P_search(:,i);
                                    %
                                    % Check constraints
                                    %
                                    if ~isempty(func_C)
                                        c_const = func_C(xtemp);
                                        m       = size(c_const,1);
                                        if sum(c_const <= 0) ~= m
                                            mask_pte(i) = false;
                                        end
                                    end
                                    %
                                    % Check if the point was already evaluated.
                                    %
                                    if mask_pte(i) && cache
                                        match = match_point(xtemp,norm_pte(i),CacheP,...
                                            CachenormP,Cache_infeasible,Cache_normInf,tol_match);
                                        if match
                                            mask_pte(i) = false;
                                        elseif i < size(P_search,2)
                                            %
                                            % Verify if next points would match this
                                            % one in the cache, if they were to be
                                            % stored one by one (every point should be
                                            % at least 'tol_match' from eachother).
                                            %
                                            index = find(abs(norm_pte(i+1:end) - norm_pte(i)) <= tol_match) + i;
                                            if ~isempty(index)
                                                auxP = P_search(:,index);
                                                index2 = max(abs(auxP-xtemp),[],1) <= tol_match;
                                                if any(index2)
                                                    index = index(index2);
                                                    mask_pte(index) = false;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            if ~all(mask_pte)
                                P_search = P_search(:,mask_pte);
                                norm_pte = norm_pte(mask_pte);
                                centre_origin = centre_origin(mask_pte);
                            end
                            if ~isempty(P_search)
                                numPointsEval = size(P_search,2);
                                if stop_feval && (func_eval + numPointsEval > max_fevals)
                                    numPointsEval = max_fevals - func_eval;
                                    halt = 1;
                                end
                                ftemp_pte = zeros(size(Flist,1), numPointsEval);
                                %
                                % Evaluate the points and store the corresponding values.
                                %
                                if parallel_version && numPointsEval > 1
                                    parfor (i = 1:numPointsEval,numPointsEval)
                                        ftemp_pte(:,i) = func_F(P_search(:,i));
                                    end
                                else
                                    for i = 1:numPointsEval
                                        ftemp_pte(:,i) = func_F(P_search(:,i));
                                    end
                                end
                                func_eval = func_eval + numPointsEval;
                                func_iter = func_iter + numPointsEval;
                                if parallel_version
                                    func_parCycles = func_parCycles + ceil(numPointsEval/numWorkers);
                                end
                                % Check if objective function values are finite
                                mask_pte = sum(isfinite(ftemp_pte),1) == size(ftemp_pte,1);
                                if halt || ~all(mask_pte)
                                    % Fill the infeasible points' cache
                                    if cache
                                        Cache_infeasible = [Cache_infeasible P_search(:,~mask_pte)];
                                        Cache_normInf = [Cache_normInf norm_pte(~mask_pte)];
                                    end
                                    % Crop results
                                    P_search = P_search(:,mask_pte);
                                    centre_origin = centre_origin(mask_pte);
                                    ftemp_pte = ftemp_pte(:,mask_pte);
                                    norm_pte = norm_pte(mask_pte);
                                    numPointsEval = size(P_search,2);
                                end
                                if ~isempty(P_search)
                                    if cache
                                        CacheP = [CacheP P_search];
                                        CachenormP = [CachenormP norm_pte];
                                        CacheF = [CacheF ftemp_pte];
                                    end
                                    for i = 1:numPointsEval
                                        [pdom,index_ndom] = paretodominance(ftemp_pte(:,i),Flist);
                                        if ~pdom
                                            if ~successful_centres(centre_origin(i))
                                                centre_idx = find(sum(Plist == centres(:,centre_origin(i)),1) == size(Plist,1), 1);
                                                successful_centres(centre_origin(i)) = true;
                                                if ~isempty(centre_idx)
                                                    alfa(centre_idx) = new_suc_alfa(centre_origin(i));
                                                end
                                            end
                                            Plist = [Plist(:,index_ndom) P_search(:,i)];
                                            Flist = [Flist(:,index_ndom) ftemp_pte(:,i)];
                                            %
                                            % Step size parameter is already
                                            % updated for the next iteration.
                                            %
                                            alfa  = [alfa(index_ndom) new_suc_alfa(centre_origin(i))];
                                            % Save new points as
                                            % centres to add for the poll step
                                            if pollOnSucSearch
                                                newPollCentres = [newPollCentres P_search(:,i)];
                                                newPollCentres_alfa = [newPollCentres_alfa new_suc_alfa(centre_origin(i))];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        level = level + 1;
                        if all(successful_centres)
                            success = 1;
                            search_success = 1;
                        end
                    end
                    if ~success && any(successful_centres)
                        success = 1;
                        search_success = 1;
                    end
                end
                % Crop successful and dominated centres; the remaining will
                % be selected as centres for the poll step
                centres = orig_centres;
                alfa_centres = orig_alfa_centres;
                if success && ~halt
                    poll = 0;
                    if numCentres > 1
                        % Check for success
                        orig_numSelected = size(orig_centres,2);
                        idx = find(cut_index);
                        idx = idx(successful_centres);
                        if length(idx) ~= orig_numSelected
                            cut_index = true(orig_numSelected,1);
                            cut_index(idx) = false;
                            centres = orig_centres(:,cut_index);
                            alfa_centres = orig_alfa_centres(cut_index);
                            numSelectedCentres = size(centres,2);
                            cut_index = true(numSelectedCentres,1);
                            % Check for dominance
                            for i = 1:numSelectedCentres
                                centre_idx = find(sum(Plist == centres(:,i),1) == size(Plist,1), 1);
                                if isempty(centre_idx)
                                    cut_index(i) = false;
                                end
                            end
                            if ~all(cut_index)
                                centres = centres(:,cut_index);
                                alfa_centres = alfa_centres(cut_index);
                            end
                            % Add the new points to our list of centres
                            % for the poll step. Now it contains the new
                            % points found in the search step + all centres
                            % that were unsuccessful at producing new
                            % points in the search step
                            if pollOnSucSearch
                                centres = [newPollCentres centres];
                                alfa_centres = [newPollCentres_alfa alfa_centres];
                                numNewCentres = size(newPollCentres,2);
                            end
                            if ~isempty(centres)
                                poll = 1;
                            end
                        end
                    end
                end
            end
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     Poll Step.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if poll && ~halt
                new_suc_alfa = alfa_centres * gamma_par;
                %
                % Generate the positive basis.
                %
                if ~dir_dense
                    D = [eye(n) -eye(n)];
                else
                    v     = 2*rand(n,1)-1;
                    [Q,R] = qr(v);
                    if ( R(1) > 1 )
                        D = Q * [ eye(n) -eye(n) ];
                    else
                        D = Q * [ -eye(n) eye(n) ];
                    end
                end
                max_D = max(sqrt(sum(D.^2)));
                %
                %        Poll using the positive basis.
                %
                numSelectedCentres = size(centres,2);
                pointsToEval = zeros(n,2*n*numSelectedCentres);
                for j = 1:numSelectedCentres
                    pointsToEval(:,(j-1)*2*n+1:j*2*n) = centres(:,j) + alfa_centres(j) * D;
                end
                numPointsEval = size(pointsToEval,2);
                mask_pte = true(1,numPointsEval);
                norm_pte = vecnorm(pointsToEval,1,1);
                % Variable that tracks the centre that originated a
                % given new point
                centre_origin = repelem(1:numSelectedCentres,2*n);
                % Variable to track success on each centre
                successful_centres = false(1,numSelectedCentres);
                for i = 1:numPointsEval
                    if mask_pte(i)
                        xtemp = pointsToEval(:,i);
                        %
                        % Check feasibility.
                        %
                        bound = [];
                        if ~isempty(ubound)
                            bound = xtemp - ubound;
                        end
                        if ~isempty(lbound)
                            bound = [bound ; - xtemp + lbound];
                        end
                        if ~isempty(bound)
                            m  = size(bound,1);
                            if sum(bound <= 0) ~= m
                                mask_pte(i) = false;
                            end
                        end
                        %
                        % Check constraints
                        %
                        if mask_pte(i) && ~isempty(func_C)
                            c_const = func_C(xtemp);
                            m       = size(c_const,1);
                            if sum(c_const <= 0) ~= m
                                mask_pte(i) = false;
                            end
                        end
                        %
                        % Check if the point was already evaluated.
                        %
                        if mask_pte(i) && cache
                            match = match_point(xtemp,norm_pte(i),CacheP,...
                                CachenormP,Cache_infeasible,Cache_normInf,tol_match);
                            if match
                                mask_pte(i) = false;
                            elseif i < numPointsEval
                                %
                                % Verify if next points would match this
                                % one in the cache, if they were to be
                                % stored one by one (every point should be
                                % at least 'tol_match' from eachother).
                                %
                                index = find(abs(norm_pte(i+1:end) - norm_pte(i)) <= tol_match) + i;
                                if ~isempty(index)
                                    auxP = pointsToEval(:,index);
                                    index2 = max(abs(auxP-xtemp),[],1) <= tol_match;
                                    if any(index2)
                                        index = index(index2);
                                        mask_pte(index) = false;
                                    end
                                end
                            end
                        end
                    end
                end
                if ~all(mask_pte)
                    pointsToEval = pointsToEval(:,mask_pte);
                    centre_origin = centre_origin(mask_pte);
                    norm_pte = norm_pte(mask_pte);
                end
                if ~isempty(pointsToEval)
                    numPointsEval = size(pointsToEval,2);
                    if stop_feval && (func_eval + numPointsEval > max_fevals)
                        numPointsEval = max_fevals - func_eval;
                        halt = 1;
                    end
                    ftemp_pte = zeros(size(Flist,1), numPointsEval);
                    %
                    % Evaluate the points and store the corresponding values.
                    %
                    if parallel_version && numPointsEval > 1
                        parfor (i = 1:numPointsEval, numPointsEval)
                            ftemp_pte(:,i) = func_F(pointsToEval(:,i));
                        end
                    else
                        for i = 1:numPointsEval
                            ftemp_pte(:,i) = func_F(pointsToEval(:,i));
                        end
                    end
                    func_eval = func_eval + numPointsEval;
                    func_iter = func_iter + numPointsEval;
                    if parallel_version
                        func_parCycles = func_parCycles + ceil(numPointsEval/numWorkers);
                    end
                    % Check if objective function values are finite
                    mask_pte = sum(isfinite(ftemp_pte),1) == size(ftemp_pte,1);
                    if halt || ~all(mask_pte)
                        % Fill the infeasible points' cache
                        if cache
                            Cache_infeasible = [Cache_infeasible pointsToEval(:,~mask_pte)];
                            Cache_normInf = [Cache_normInf norm_pte(~mask_pte)];
                        end
                        % Crop results
                        pointsToEval = pointsToEval(:,mask_pte);
                        ftemp_pte = ftemp_pte(:,mask_pte);
                        centre_origin = centre_origin(mask_pte);
                        norm_pte = norm_pte(mask_pte);
                        numPointsEval = size(pointsToEval,2);
                    end
                    if ~isempty(pointsToEval)
                        if cache
                            CacheP = [CacheP pointsToEval];
                            CachenormP = [CachenormP norm_pte];
                            CacheF = [CacheF ftemp_pte];
                        end
                        for i = 1:numPointsEval
                            [pdom,index_ndom] = paretodominance(ftemp_pte(:,i),Flist);
                            if ~pdom
                                if ~successful_centres(centre_origin(i))
                                    centre_idx = find(sum(Plist == centres(:,centre_origin(i)),1) == size(Plist,1), 1);
                                    successful_centres(centre_origin(i)) = true;
                                    if ~isempty(centre_idx)
                                        alfa(centre_idx) = new_suc_alfa(centre_origin(i));
                                    end
                                end
                                success = 1;
                                Plist = [Plist(:,index_ndom), pointsToEval(:,i)];
                                Flist = [Flist(:,index_ndom), ftemp_pte(:,i)];
                                %
                                % Step size parameter is already
                                % updated for the next iteration.
                                %
                                alfa  = [alfa(index_ndom) new_suc_alfa(centre_origin(i))];
                            end
                        end
                    end
                end
            end
            if ~halt
                numOrigCentres = size(orig_centres,2);
                index_poll_centres = zeros(1,numOrigCentres + numNewCentres);
                if search_success && poll
                    idxFromSearch = zeros(1,numOrigCentres);
                    idxFromPoll = zeros(1,numOrigCentres + numNewCentres);
                end
                %
                % Verify which original centers are maintained in the end;
                % the first cycle is for the new points generated in the
                % search step (if any), the second cycle is for the centres
                % selected at the beggining.
                %
                for k = 1:numNewCentres
                    centre_idx = find(sum(Plist == centres(:,k),1) == size(Plist,1), 1);
                    if ~isempty(centre_idx)
                        index_poll_centres(k) = centre_idx;
                        idxFromPoll(k) = centre_idx;
                    end
                end
                %
                % If the search step was successful and the poll step was also
                % required, place the successful centers' indexes first
                % on the block to rotate [centres_search | centres_poll]
                %
                for k = 1:numOrigCentres
                    centre_idx = find(sum(Plist == orig_centres(:,k),1) == size(Plist,1), 1);
                    if ~isempty(centre_idx)
                        index_poll_centres(k + numNewCentres) = centre_idx;
                        if search_success && poll
                            poll_idx = find(sum(centres == orig_centres(:,k),1) == size(Plist,1), 1);
                            if isempty(poll_idx)
                                idxFromSearch(k) = centre_idx;
                            else
                                idxFromPoll(k + numNewCentres) = centre_idx;
                            end
                        end
                    end
                end
                index_poll_centres = index_poll_centres(logical(index_poll_centres));
                if search_success && poll
                    idxFromSearch = idxFromSearch(logical(idxFromSearch));
                    idxFromPoll = idxFromPoll(logical(idxFromPoll));
                end
                %
                % Update the counter for successful iterations
                %
                % Update the step size on the original centres that were
                % not successful (successful iteration)
                %
                % Update the step size on all original centres
                % (unsuccessful iteration)
                %
                if success
                    iter_suc = iter_suc + 1;
                    % Retrieve unsuccessful points and reduce step size
                    if poll && ~isempty(index_poll_centres)
                        for i = 1:numSelectedCentres
                            if ~successful_centres(i)
                                centre_idx = find(sum(Plist == centres(:,i),1) == size(Plist,1), 1);
                                if ~isempty(centre_idx)
                                    alfa(centre_idx) = alfa(centre_idx) * beta_par;
                                end
                            end
                        end
                    end
                elseif ~isempty(index_poll_centres)
                    alfa(index_poll_centres) = alfa(index_poll_centres) * beta_par;
                end
                %
                % Rotate maintained centres to the end of the list
                %
                if ~isempty(index_poll_centres)
                    nPlist = size(Plist,2);
                    tmp_allidx = 1:nPlist;
                    tmp_allidx(index_poll_centres) = 0;
                    first_idx = tmp_allidx(logical(tmp_allidx));
                    if search_success && poll
                        second_idx = [idxFromSearch idxFromPoll];
                    else
                        second_idx = index_poll_centres;
                    end
                    Plist  = [Plist(:,first_idx) Plist(:,second_idx)];
                    Flist  = [Flist(:,first_idx) Flist(:,second_idx)];
                    alfa   = [alfa(first_idx) alfa(second_idx)];
                end
            elseif success
                iter_suc = iter_suc + 1;
            end
        else
            %
            % If the first point of the list is below the alfa stopping
            % criterion, rotate it to the end of the list.
            %
            nPlist = size(Plist,2);
            Plist  = [Plist(:,2:nPlist),Plist(:,1)];
            Flist  = [Flist(:,2:nPlist),Flist(:,1)];
            alfa   = [alfa(2:nPlist) alfa(1)];
            move   = 1;
        end
        if ~move
            %
            % Check if the stopping criteria are satisfied.
            %
            if stop_alfa && (sum(alfa >= tol_stop)==0)
                halt = 1;
            end
            if stop_feval && (func_eval >= max_fevals)
                halt = 1;
            end
            if parallel_version && stop_fparcycles && (func_parCycles >= max_fparCycles)
                halt = 1;
            end
            iter = iter + 1;
            %
            % Print the iteration report.
            %
            print_format = ['| %5d |       %2d       |    %2d   |   %5d   | %+13.8e  | %+13.8e  |\n'];
            fprintf(print_format, iter, search_success, success, size(Plist,2), min(alfa), max(alfa));
            fprintf(fresult,print_format, iter, search_success, success, size(Plist,2), min(alfa), max(alfa));
            %
            % Save last iteration values for unexpected error recovery.
            %
            last_Plist  = Plist;
            last_Flist  = Flist;
            last_alfa   = alfa;
            last_iter_suc = iter_suc;
            last_func_eval = func_eval;
            if parallel_version
                last_func_parCycles = func_parCycles;
            end
            last_max_D = max_D;
            if cache
                last_CacheP     = CacheP;
                last_CachenormP = CachenormP;
                last_CacheF     = CacheF;
                last_Cache_infeasible = Cache_infeasible;
                last_Cache_normInf = Cache_normInf;
            end
            last_gamma_cycle = gamma_cycle;
        end
    end
    time = etime(clock,time);
    if list == 4
        time = time + prev_time;
    end
    prev_time = time;
    time_set = true;
    list_size = size(Plist,2);
    %
    % Print final report in screen.
    %
    fprintf('\n Final Report: \n\n');
    print_format = 'Elapsed Time = %10.3e \n\n';
    fprintf(print_format,time);
    if parallel_version
        fprintf('| #iter | #isuc | list size | #fevals | #fparCycles |\n');
        print_format = ['| %5d | %5d |   %5d   |  %5d  |   %5d     |\n\n'];
        fprintf(print_format, iter, iter_suc,list_size,func_eval,func_parCycles);
    else
        fprintf('| #iter | #isuc | list size | #fevals |\n');
        print_format = ['| %5d | %5d |   %5d   |  %5d  |\n\n'];
        fprintf(print_format, iter, iter_suc,list_size,func_eval);
    end
    %
    % Print final report in the results file.
    %
    fprintf(fresult,'\n Final Report: \n\n');
    print_format = 'Elapsed Time = %10.3e \n\n';
    fprintf(fresult,print_format,time);
    if parallel_version
        fprintf(fresult,'| #iter | #isuc | list size | #fevals | #fparCycles |\n');
        print_format = ['| %5d | %5d |   %5d   |  %5d  |   %5d     |\n\n'];
        fprintf(fresult,print_format, iter, iter_suc,list_size,func_eval,func_parCycles);
    else
        fprintf(fresult,'| #iter | #isuc | list size | #fevals |\n');
        print_format = ['| %5d | %5d |   %5d   |  %5d  |\n\n'];
        fprintf(fresult,print_format, iter, iter_suc,list_size,func_eval);
    end
    fclose(fresult);
    %
    % Save last iteration matrix results to file
    %
    cd(foldername);
    variableFile = "BoostDMS_lastiteration_" + problemName + ".mat";
    if cache
        save(variableFile,"Plist","Flist","alfa","CacheP","CachenormP","CacheF","prev_time",...
            "iter","iter_suc","func_eval","func_parCycles","max_D","lbound","ubound",...
            "gamma_cycle","Cache_infeasible","Cache_normInf");
    else
        save(variableFile,"Plist","Flist","alfa","prev_time",...
            "iter","iter_suc","func_eval","func_parCycles","max_D","lbound","ubound",...
            "gamma_cycle");
    end
    cd(fullfile("../"));
    %
    %
    % Print a plot of the computed Pareto front.
    %
    plotname = "BoostDMS_plot_" + problemName + ".jpg";
    figname = "BoostDMS_plot_matlab_" + problemName + ".fig";
    if size(Flist,1) == 2
        cd(foldername);
        h = figure;
        plot(Flist(1,:),Flist(2,:),'LineStyle','None','Marker','o');
        xlabel('f1');
        ylabel('f2');
        savefig(figname);
        print(h,'-djpeg', plotname);
        cd(fullfile("../"));
    elseif size(Flist,1) == 3
        cd(foldername);
        h = figure;
        scatter3(Flist(1,:),Flist(2,:),Flist(3,:),'Marker','o');
        xlabel('f1');
        ylabel('f2');
        zlabel('f3');
        savefig(figname);
        print(h,'-djpeg', plotname);
        cd(fullfile("../"));
    end
    %
catch exception
    if iter > 0 && iter > start_iter
        Plist = last_Plist;
        Flist = last_Flist;
        alfa = last_alfa;
        if ~time_set
            time = etime(clock,time);
            if list == 4
                time = time + prev_time;
            end
        end
        prev_time = time;
        iter_suc = last_iter_suc;
        func_eval = last_func_eval;
        if parallel_version
            func_parCycles = last_func_parCycles;
        end
        max_D = last_max_D;
        if cache
            CacheP = last_CacheP;
            CachenormP = last_CachenormP;
            CacheF = last_CacheF;
            Cache_infeasible = last_Cache_infeasible;
            Cache_normInf = last_Cache_normInf;
        end
        gamma_cycle = last_gamma_cycle;
        %
        % Save last iteration matrix results to file
        %
        cd(foldername);
        variableFile = "BoostDMS_lastiteration_" + problemName + ".mat";
        if cache
            save(variableFile,"Plist","Flist","alfa","CacheP","CachenormP","CacheF","prev_time",...
                "iter","iter_suc","func_eval","func_parCycles","max_D","lbound","ubound",...
                "gamma_cycle","Cache_infeasible","Cache_normInf");
        else
            save(variableFile,"Plist","Flist","alfa","prev_time",...
                "iter","iter_suc","func_eval","func_parCycles","max_D","lbound","ubound",...
                "gamma_cycle");
        end
        cd(fullfile("../"));
    end
    rethrow(exception);
end
%
% End BoostDMS.