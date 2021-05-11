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
grad_C = zeros(8,3);
% -----------------------------------------------------------------------------
%  Block to be user modified.
% -----------------------------------------------------------------------------
r_b = x(1); r_p = x(2); d_b = x(3); d_p = x(4); d = x(5); h = x(6); phi = x(7); beta = x(8); 
% C(1) = (r_p*cos(-d_p) - r_b*cos(-d_b) - h*cos(-beta))^2 + (r_p*sin(-d_p) - r_b*sin(-d_b) - h*sin(-beta))^2 - d^2;
% C(2) = r_b*sin(-d_b) + h*sin(-beta) - r_b*sin(d_b)-h*sin(beta);
% C(3) = r_b*sin(4/3*pi + d_b) + h*sin(4/3*pi + beta) + r_b*sin(+d_b) + h*sin(+beta);

%diff((r_p*cos(-d_p) - r_b*cos(-d_b) - h*cos(-beta))^2 + (r_p*sin(-d_p) - r_b*sin(-d_b) - h*sin(-beta))^2 - d^2, var)
grad_C(1,1) = 2*cos(d_b)*(h*cos(beta) + r_b*cos(d_b) - r_p*cos(d_p)) + 2*sin(d_b)*(h*sin(beta) + r_b*sin(d_b) - r_p*sin(d_p));
grad_C(2,1) = - 2*cos(d_p)*(h*cos(beta) + r_b*cos(d_b) - r_p*cos(d_p)) - 2*sin(d_p)*(h*sin(beta) + r_b*sin(d_b) - r_p*sin(d_p));
grad_C(3,1) = 2*r_b*cos(d_b)*(h*sin(beta) + r_b*sin(d_b) - r_p*sin(d_p)) - 2*r_b*sin(d_b)*(h*cos(beta) + r_b*cos(d_b) - r_p*cos(d_p));
grad_C(4,1) = 2*r_p*sin(d_p)*(h*cos(beta) + r_b*cos(d_b) - r_p*cos(d_p)) - 2*r_p*cos(d_p)*(h*sin(beta) + r_b*sin(d_b) - r_p*sin(d_p));
grad_C(5,1) = -2*d;
grad_C(6,1) = 2*cos(beta)*(h*cos(beta) + r_b*cos(d_b) - r_p*cos(d_p)) + 2*sin(beta)*(h*sin(beta) + r_b*sin(d_b) - r_p*sin(d_p));
grad_C(7,1) = 0;
grad_C(8,1) = 2*h*cos(beta)*(h*sin(beta) + r_b*sin(d_b) - r_p*sin(d_p)) - 2*h*sin(beta)*(h*cos(beta) + r_b*cos(d_b) - r_p*cos(d_p));

grad_C(1,3) = -2*sin(d_b);
grad_C(2,2) = 0;
grad_C(3,2) = -2*r_b*cos(d_b);
grad_C(4,2) = 0;
grad_C(5,2) = 0;
grad_C(6,2) = -2*sin(beta);
grad_C(7,2) = 0;
grad_C(8,2) = -2*h*cos(beta);

grad_C(1,3) = sin(4*pi/3 + d_b) + sin(d_b);
grad_C(2,3) = 0;
grad_C(3,3) = cos(4/3*pi + d_b) + cos(d_b);
grad_C(4,2) = 0;
grad_C(5,2) = 0;
grad_C(6,2) = sin(4/3*pi + d_b) + sin(beta);
grad_C(7,2) = 0;
grad_C(8,2) = h*(cos(4/3*pi + d_b) + h*cos(beta));
%
% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
% End of grad_const.