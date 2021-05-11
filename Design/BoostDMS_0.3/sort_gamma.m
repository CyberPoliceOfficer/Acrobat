function [Psort,Fsort,alfasort,gamma_cycle] = sort_gamma(P,F,alfa,stop_alfa,tol_stop,spread_option,gamma_cycle)
%
% Purpose:
%
%    Function sort_gamma sorts a list of points according to the largest
%    gap between two consecutive points in the list considered.
%
% Input:
%
%         P (List of points.)
%
%         F (Function values corresponding to the points in the list.)
%
%         alfa (Step size parameters corresponding to the points in the list.)
%
%         stop_alfa (0-1 variable, which indicates if there is a stopping criterion
%                    based on the step size.)
%
%         tol_stop (Lowest value allowed for the step size parameter.)
%
%         spread_option (0-1 variable, which indicates the type of ordering
%                      strategy to consider.)
%
%         gamma_cycle (Variable that defines the first component to be
%           selected, in case of cyclic component selection.)
%
% Output:
%
%         Psort (Sorted list of points.)
%
%         Fsort (Function values corresponding to the points in the sorted list.)
%
%         alfasort (Step size parameters corresponding to the points in the
%         sorted list.)
%
%         gamma_cycle (Variable that defines the first component to be
%           selected, in case of cyclic component selection.)
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
Psort    = P;
Fsort    = F;
alfasort = alfa;
if stop_alfa
    index1 = alfa>=tol_stop;
else
    index1 = true;
end
if ~all(index1)
    Paux   = P(:,index1);
    Faux   = F(:,index1);
    alfaux = alfa(index1);
    nPaux  = size(Paux,2);
else
    Paux   = P;
    Faux   = F;
    alfaux = alfa;
    nPaux  = size(Paux,2);
end
if nPaux >= 2
    [alfaux,index_alfa] = sort(alfaux,'descend');
    Paux           = Paux(:,index_alfa);
    Faux           = Faux(:,index_alfa);
    %       if spread_option == 0 || spread_option == 1
    nobj = size(Faux,1);
    if spread_option == 1
        Mdist = -inf*ones(2,nPaux);
        for_range = gamma_cycle;
        gamma_cycle = mod(gamma_cycle,nobj) + 1;
    else
        Mdist = -inf*ones(2*nobj,nPaux);
        for_range = 1:nobj;
    end
    for j = for_range
        [FProj,index] = sort(Faux(j,:));
        if spread_option == 1
            k = 1;
        else
            k = j;
        end
        Mdist(2*k,index(1)) = FProj(2)-FProj(1);
        for i = 2:nPaux-1
            Mdist(2*k-1,index(i)) = FProj(i)-FProj(i-1);
            Mdist(2*k,index(i))   = FProj(i+1)-FProj(i);
        end
        Mdist(2*k-1,index(nPaux)) = FProj(nPaux)-FProj(nPaux-1);
    end
    [~,index] = sort(max(Mdist),'descend');
    %       end
    Paux   = Paux(:,index);
    alfaux = alfaux(index);
    if stop_alfa
        index2 = alfa<tol_stop;
    else
        index2 = 0;
    end
    if ~any(index2)
        Psort    = Paux;
        Fsort    = F(:,index_alfa(index));
        alfasort = alfaux;
    else
        Psort    = [Paux,P(:,index2)];
        Fsort    = F(:,index1);
        Fsort    = [Fsort(:,index_alfa(index)),F(:,index2)];
        alfasort = [alfaux,alfa(index2)];
    end
end
%
% End of sort_gamma.