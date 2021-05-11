function  [quad,P_search,cut_index,g_old,H_old] = individual_models(tol_match,Delta,...
    P,F,centres)
%
% Purpose:
%
%        Function individual_models computes a quadratic polynomial model
%        for each component of the objective function and provides a list
%        of points corresponding to the individual minimization of each model
%        inside the trust region considered. Several points may be
%        considered as model centres.
%
% Input:
%
%         tol_match (tolerance used in point comparisons, when considering
%                   a cache.)
%
%         Delta (radius of the trust region, for every provided centre.)
%
%         P (list of points to be used in the computation of the models,
%           for every provided centre.)
%
%         F (function values corresponding to the points in the list,
%           for every provided centre.)
%
%         centres (set of points to be considered as model centres.)
%
% Output:
%
%         quad      (0-1 variable: 1 if a quadratic model is computed;
%                       0 otherwise.)
%
%         P_search  (list of candidate points to be evaluated.)
%
%         cut_index (logical vector recording the centres that generated
%           valid models.)
%
%         g_old     (structure storing the gradients of the models for all
%         centres.)
%
%         H_old     (structure storing the Hessians of the models for all
%         centres.)
%
% Functions called: quad_Frob
%
% BoostDMS Version 0.2.
% Copyright (C) 2020 C. P. Brás, A. L. Custódio, V. M. Duarte, P. D. Medeiros, S. Tavares
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
P_search = [];
g_old = [];
H_old = [];
quad = 0;
if iscell(P)
    numCentres = length(P);
else
    numCentres = 1;
end
cut_index = true(numCentres,1);
for i = 1:numCentres
    if Delta(i) <= tol_match
        cut_index(i) = false;
    end
end
%
% Compute the quadratic models.
%
if any(cut_index)
    %
    % Check if there are enough points in P.
    %
    numCentres = sum(cut_index);
    n = size(centres,1);
    if iscell(P)
        for i = 1:numCentres
            tmp_ind = find(cut_index,i);
            index = tmp_ind(end);
            if size(P{index},2) <= n+1
                cut_index(index) = false;
            end
        end
        if ~all(cut_index)
            P = P(cut_index);
            centres = centres(:,cut_index);
            Delta = Delta(cut_index);
            numCentres = length(P);
        end
        if ~isempty(P)
            quad = 1;
        end
    elseif size(P,2) > n+1
        quad = 1;
    end
    %
    if quad
        if iscell(F)
            models_per_centre = size(F{1},1);
        else
            models_per_centre = size(F,1);
        end
        P_search = zeros(n,models_per_centre,numCentres);
        g_old = zeros(n,models_per_centre,numCentres);
        H_old = zeros(n,n,models_per_centre,numCentres);
        for i = 1:numCentres
            for j = 1:models_per_centre
                %
                % Build the model.
                %
                if iscell(F)
                    [Haux,gaux]  = quad_Frob(P{i},F{i}(j,:));
                else
                    [Haux,gaux]  = quad_Frob(P,F(j,:));
                end
                H_old(:,:,j,i) = Haux;
                g_old(:,j,i) = gaux;
                %
                % Minimize the model in the trust region.
                %
                if norm(gaux) <= eps
                    [U,e]      = eig(Haux);
                    [~,indexe] = sort(diag(e));
                    P_search(:,j,i) = centres(:,i) + Delta(i)* U(:,indexe(1))/norm(U(:,indexe(1)));
                else
                    P_search(:,j,i) = centres(:,i) + trust(gaux,Haux,Delta(i));
                end
            end
        end
        if numCentres > 1
            P_search = reshape(P_search, n, models_per_centre * numCentres);
        end
    end
end
%
% End of individual_models.