function [f,g] = obj_model(x)
%
% Purpose:
%
%    Function obj_model computes the objective function value and the gradient
%    vector at point x, corresponding to the following problem to be solved 
%    by fmincon:
%
%    min  zeta
%    s.t. m_j(x_center + x) <= zeta, j in a subset of the components of a 
%                                    multiobjective function 
%         ||x||  <= Delta
%         lbound <= x_center + x <= ubound   
%
%    where m_j(x_center + x) is a quadratic polynomial model defined around
%    x_center for the j component of a multiojective function.
%       
% Input:
%
%    x  (point where the objective function is evaluated. The last position  
%       is the value of variable zeta.)
%          
% Output:
%
%    f  (function value at x.)
%    g  (gradient vector at x.)
%
% BoostDMS Version 0.2.
% Copyright (C) 2020 C. P. Brás, A. L. Custódio, V. M. Duarte, P. D. Medeiros, S. Tavares
% http://ferrari.dmat.fct.unl.pt/personal/alcustodio/BoostDFO.htm
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA,
% or see <https://www.gnu.org/licenses/>.
%
%
n = length(x);
f = x(n);
if nargout > 1
    g = [zeros(n-1,1);1];
end
%
% End of obj_model.



