function [grad_C] = grad_C(x);
%
% Purpose:
%
%    Function grad_C is an user provided function which 
%    computes the gradients of the constraints functions, other
%    than bounds, c_j(x) <= 0, j = 1,...,p, at a point provided 
%    by the optimizer. Bounds on the variables, should be
%    indicated in the BoostDMS command line.
%
%
% Input:  
%
%         x (Point given by the optimizer.)
%
% Output: 
%
%         grad_C (Matrix storing columnwise the gradients of the functions
%                 c_j defining the constraints at the given point.)
%
%
% SID-PSM Version 1.3
% Copyright (C) 2009 A. L. Custodio and L. N. Vicente.
% http://www.mat.uc.pt/sid-psm
%
% BoostDMS Version 0.1
% Copyright (C) 2019 C. P. Br�s and A. L. Cust�dio.
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
grad_C = zeros(24,20);
% -----------------------------------------------------------------------------
%  Block to be user modified.
% -----------------------------------------------------------------------------
j = 1;
for i = 2:24
    if (rem(i,6) ~= 1)
        grad_C(i,j) = 1;
        grad_C(i-1,j) = -1;
        j = j+1;
    end
end
%
% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
% End of grad_const.