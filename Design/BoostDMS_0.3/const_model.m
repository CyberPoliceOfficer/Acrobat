function [c,ceq,gradc,gradceq] = const_model(x,v,Fx,gi,Hi,Delta,func_C,grad_C)
%
% Purpose:
%
%    Function const_model computes the nonlinear constraints and the 
%    gradients at point x, corresponding to the following problem to be
%    solved by fmincon:
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
%    x     (point where the constraints functions are evaluated. The last position  
%          is the value of variable zeta.)
%    v     (indexes of components of the multiobjective functions to be
%          considered.)
% 	 Fx    (values of all components of the objective function at x_center.)
%    gi    (structure with all objective gradients of the quadratic models.)
%    Hi    (structure with all objective Hessians of the quadratic models.)
%    Delta (radius of the trust region.)
%    func_C (Name of the file for computing other problem constraints.)
%    grad_C (Name of the file for computing gradients of problem constraints,
%           other than bounds.)
%
% Output:
%
%    c, ceq, gradc, gradceq (nonlinear inequalities or equalities and
%                           the respective gradients evaluated at x).   
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
n    = size(x,1);
x1   = x(1:n-1);
nobj = size(v,2);
%
ceq  = [];
c    = zeros(nobj+1,1);
for i = 1:nobj
    c(i) = Fx(v(i))+ x1'*gi(:,v(i))+0.5*x1'*Hi(:,:,v(i))*x1 - x(n);
end
c(nobj+1) = x1'*x1 - Delta^2;
%
if ~isempty(func_C)
    c_Omega = func_C(x1);
    const   = length(c_Omega);
    c       = cat(1,c,c_Omega);
end
%
if nargout > 2
   gradceq = [];
   gradc   = zeros(n,nobj+1);
   for i = 1:nobj
      gradc(1:n-1,i) = gi(1:n-1,v(i))+Hi(1:n-1,1:n-1,v(i))*x1;
   end
   gradc(n,:)          = [-ones(1,nobj),0];
   gradc(1:n-1,nobj+1) = 2*x1;
   if ~isempty(func_C)
      grad_Omega      = grad_C(x1);
      grad_Omega(n,:) = zeros(1,const);
      gradc           = cat(2,gradc,grad_Omega);
   end
end
%
% End of const_model.
    
    