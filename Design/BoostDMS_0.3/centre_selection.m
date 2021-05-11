function [centres,new_alfa] = centre_selection(Plist,Flist,alfa,stop_alfa,tol_stop,numCentres)
%
% Purpose:
%
%    Function centre_selection selects all model centers according
%       to the selection strategy chosen.
%
% Input:
%
%         Plist (List storing columnwise the nondominated points.)
%
%         Flist (Function values corresponding to the points in 'Plist'.)
%
%         alfa (alfa values corresponding to the points in 'Plist'.)
%
%         stop_alfa (0-1 variable, which indicates if there is a stopping
%           criterion based on the step size.)
%
%         tol_stop (Lowest value allowed for the step size parameter.)
%
%         numCentres (Number of centres to select.)
%
% Output:
%
%         centres (List storing columnwise the points chosen as model centres.)
%
%         new_alfa (List of alfa values that correspond to each of
%            the centres found, in case the iteration is found successful).
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
centresFromSpread = 2;
if stop_alfa
    index1 = alfa >= tol_stop;
    if ~all(index1)
        Plist = Plist(:,index1);
        Flist = Flist(:,index1);
        alfa  = alfa(index1);
    end
end
if size(Plist,2) <= numCentres
    centres = Plist;
    new_alfa = alfa;
else
    cut_index = true(1,size(Plist,2));
    nobj = size(Flist,1);
    indexes = 1:numCentres;
    if centresFromSpread
        cut_index(1:centresFromSpread) = false;
    end
    selected = centresFromSpread;
    for i = 1:nobj
        [m,index] = min(Flist(i,cut_index));
        m_values_idx = find(Flist(i,cut_index) == m);
        % In case of draw, select the biggest alfa
        if length(m_values_idx) == 1
            tmpindex = find(cut_index,index);
            index = tmpindex(end);
        else
            tmpalfas = zeros(length(m_values_idx),1);
            tmpindexes = zeros(length(m_values_idx),1);
            for k = 1:length(m_values_idx)
                tmp_idx = find(cut_index,m_values_idx(k));
                tmpindexes(k) = tmp_idx(end);
                tmpalfas(k) = alfa(tmp_idx(end));
            end
            [~,index] = max(tmpalfas);
            index = tmpindexes(index);
        end
        selected = selected + 1;
        indexes(selected) = index;
        cut_index(index) = false;
    end
    centres = Plist(:,indexes);
    new_alfa = alfa(indexes);
end
%
% End of centre_selection.