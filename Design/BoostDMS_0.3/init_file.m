%
%    This file contains all the initialization variables needed to 
%       solve the problem with BoostDMS. These include the initial 
%       point(s) and bounds. Bounds are always required, x_initial
%       may be empty, a point or a list of points, depending on the 
%       'list' option selected on 'parameters_BoostDMS'.
%     
%
% BoostDMS Version 0.2
% Copyright (C) 2020 C. P. Br�s, A. L. Cust�dio, V. M. Duarte, P. D. Medeiros, S. Tavares
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

%r_b = x(1); r_p = x(2); d_b = x(3); d_p = x(4); d = x(5); h = x(6); phi = x(7); beta = x(8);
x_initial = [];
lowerbound = [0.01; 0.01; 0; 0; 0.01; 0.01; -pi/3; 0];
upperbound = [0.06; 0.06; pi/3; pi/3; 0.2; 0.2; pi/3; pi];

% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%