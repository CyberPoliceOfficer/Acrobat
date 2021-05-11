function match = match_point(x,xnorm,CacheP,CachenormP,Cache_infeasible,Cache_normInf,tol_match)
%
% Purpose:
%
%    Function match_point scans a list of previously evaluated points
%    to try to match a point provided by the optimizer.
%
% Input:
%
%         x (point to be checked.)
%
%         xnorm (1-norm of the point to be checked.)
%
%         CacheP (matrix of points to be used in point comparisons,
%                storing the points columnwise.)
%
%         CachenormP (vector storing the 1-norm of the points in CacheP.)
%
%         Cache_infeasible (same as CacheP, containing infeasible points 
%           from hidden constraints.)
%
%         Cache_normInf (same as CachenormP, containing norms of infeasible
%           points from hidden constraints.)
%
%         tol_match (tolerance value within which two points are
%                   considered as equal.)
%
% Output:
%
%         match (0-1 variable: 1 if x was previously evaluated; 0
%         otherwise.)
%
% DMS Version 0.2.
% Copyright (C) 2011 A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz,
% and L. N. Vicente.
% http://www.mat.uc.pt/dms
%
% BoostDMS Version 0.3.
% Copyright (C) 2020 C. P. Brás, A. L. Custódio, V. M. Duarte, P. D. Medeiros, S. Tavares
% http://ferrari.dmat.fct.unl.pt/personal/alcustodio/BoostDFO.htm
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
%
%
% Prune the points that do not satisfy a 1-norm criterion.
%
match = 0;
%
% Test if the point is in the regular cache
%
index = find(abs(CachenormP - xnorm) <= tol_match);
if ~isempty(index)
    CacheP = CacheP(:,index);
    %
    % Finish search.
    %
    index2 = find(max(abs(CacheP-x),[],1) <= tol_match,1);
    match = ~isempty(index2);
end
%
% Test if the point is in the infeasible cache
%
if ~match && ~isempty(Cache_infeasible)
    index = find(abs(Cache_normInf - xnorm) <= tol_match);
    if ~isempty(index)
        Cache_infeasible = Cache_infeasible(:,index);
        index2 = find(max(abs(Cache_infeasible-x),[],1) <= tol_match,1);
        match = ~isempty(index2);
    end
end
%
% End of match_point.