function [F] = func_F(x)
%
% Purpose:
%
%    Function func_F is an user provided function which
%    computes the value of the objective function at a
%    point provided by the optimizer.
%
% Input:  
%
%         x (Point given by the optimizer.)
%
% Output: 
%
%         F (Function value at the given point.)        
%
% DMS Version 0.2.
% BoostDMS Version 0.1.
%
% Copyright (C) 2011 A. L. Cust�dio, J. F. A. Madeira, A. I. F. Vaz, 
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
% -----------------------------------------------------------------------------
%  Block to be user modified.
% -----------------------------------------------------------------------------
F = zeros(2,1); % allocate output
r = [cos(x(7:12)).*sin(x(1:6)); sin(x(7:12)).*sin(x(1:6)); cos(x(1:6))];
u = [cos(x(19:24)).*sin(x(13:18)); sin(x(19:24)).*sin(x(13:18)); cos(x(13:18))];
A = [u; cross(r,u)];

S = svd(A);
kinv = S(6)/S(1);

if (kinv > 1e-8)
    Ainv = inv(A);
    b = Ainv(:,1:3);
    c = Ainv(:,4:6);
    F(1) = -min(1./vecnorm(b'));
    F(2) = -min(1./vecnorm(c'));
end
%-----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
% End of func_F.