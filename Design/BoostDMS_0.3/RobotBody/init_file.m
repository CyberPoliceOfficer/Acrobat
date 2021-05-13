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

x_initial = [];
lowerbound = zeros(1,24);
upperbound = ones(1,24);
upperbound(1:6) = pi;
upperbound(7:12) = 2*pi;
upperbound(13:18) = pi;
upperbound(19:24) = 2*pi;

% -----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%