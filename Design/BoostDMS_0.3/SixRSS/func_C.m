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
C = zeros(3,1);
%
% -----------------------------------------------------------------------------
%  Block to be user modified.
% -----------------------------------------------------------------------------
%     %% Compute vector bk and pk
% 
%     k = [1 2 3 4 5 6];
%     n = floor((k-1)/2);
% 
%     theta_b = n*(2/3)*pi + (-1).^k*d_b;
%     theta_p = n*(2/3)*pi + (-1).^k*d_p;
% 
%     b_k = [r_b * cos(theta_b); r_b * sin(theta_b);zeros(1,6)];
%     p_k = [r_p * cos(theta_p); r_p * sin(theta_p);zeros(1,6)];
% 
% 
%     %% Compute beta and gamma
% 
%     beta_k = n*(2/3)*pi + (-1).^(k)*beta;
%     phi_k = (-1).^(k+1)*phi;

r_b = x(1); r_p = x(2); d_b = x(3); d_p = x(4); d = x(5); h = x(6); phi = x(7); beta = x(8);
theta_b = -d_b;
theta_p = -d_p;
beta_k = -beta;   
C(1) = (r_p*cos(theta_p) - r_b*cos(theta_b) - h*cos(beta_k)).^2 + (r_p*sin(theta_p) - r_b*sin(theta_b) - h*sin(beta_k)).^2 - d^2;
C(2) = -2*(r_b*sin(d_b) + h*sin(beta));
C(3) = r_b*sin(4/3*pi + d_b) + h*sin(4/3*pi + beta) + r_b*sin(d_b) + h*sin(beta);

% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
%C = C';
%
% End of func_C.