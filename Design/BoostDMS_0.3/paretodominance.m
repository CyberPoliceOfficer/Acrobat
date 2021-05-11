function [pdom,index_ndom] = paretodominance(F,Flist)
%
% Purpose:
%
%    Function paretodominance checks if the vector F satisfies a Pareto
%    dominance criterion.
%
% Input:  
%
%         F (Vector to be checked.)
%
%         Flist (List storing columnwise the objective function values 
%                at nondominated points.)
%
% Output: 
%
%         pdom (0-1 variable: 1 if the point is dominated; 0 otherwise.)
%
%         index_ndom (0-1 vector: 1 if the corresponding list point is not
%                    dominated; 0 otherwise.)
%
% DMS Version 0.2.
% Copyright (C) 2011 A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, 
% and L. N. Vicente.
% http://www.mat.uc.pt/dms
%
% BoostDMS Version 0.2.
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
index_ndom = [];
pdom = 1;
if all(sum(F<Flist,1))
    pdom       = 0;
    index_ndom = sum(F<=Flist,1) ~= size(Flist,1);
end
%
% End of paretodominance.