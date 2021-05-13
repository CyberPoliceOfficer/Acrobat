function [C] = func_C(x)
%
% Purpose:
%
%    Function func_C is an user provided function which
%    computes the values of the constraints functions, other
%    than bounds, c_j(x) <= 0, j = 1,...,p, at a point provided 
%    by the optimizer. Bounds on the variables, should be
%    indicated in the BoostDMS command line.
%
% Input:  
%
%         x (Point given by the optimizer.)
%
% Output: 
%
%         C (Vector storing the values of the functions c_j
%           defining the constraints at the given point.)
%
% DMS Version 0.2.
% BoostDMS Version 0.1.
%
% Copyright (C) 2011 A. L. Cust???dio, J. F. A. Madeira, A. I. F. Vaz, 
% and L. N. Vicente.
%
% http://www.mat.uc.pt/dms
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
C = zeros(20,1);
%
% -----------------------------------------------------------------------------
%  Block to be user modified.
% -----------------------------------------------------------------------------

j=1;
for i = 2:24
    if (rem(i,6) ~= 1)
        C(j) = x(i) - x(i-1);
        j=j+1;
    end
end

% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
%C = C';
%
% End of func_C.