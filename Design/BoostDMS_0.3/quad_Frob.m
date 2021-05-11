function [H,g] = quad_Frob(X,F_values)
%
% Purpose:
%
%    Function quad_Frob computes a quadratic interpolation model by 
%    minimizing the Frobenius norm of the Hessian matrix. When the number
%    of points in the interpolation set exceeds the total number of points
%    required for computing a complete quadratic model, regression techniques
%    are considered.
%
% Input:  
%
%        X (matrix storing the points columnwise).
%
%        F_values (vector storing the objective function values at the
%        points).
%
% Output: 
%
%        H (matrix storing the Hessian of the model).
%
%        g (vector storing the model gradient).
%
% DMS Version 0.2.
% Copyright (C) 2011 A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, 
% and L. N. Vicente.
% http://www.mat.uc.pt/dms
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
% Initialization step.
%
tol_svd = eps;   % Minimum value accepted for a singular value.
%
[n,m] = size(X);
H     = zeros(n);
%
Y        = X;
Y_values = F_values;
Y        = Y - diag(X(:,1))*ones(n,m);
%
% Compute a quadratic model by minimizing the Frobenius norm of the Hessian.
%
if (m <= (n+1)*(n+2)/2)
    b = [Y_values zeros(1,n+1)]';
    A = ((Y'*Y).^2)/2;
    W = [A ones(m,1) Y';[ones(1,m);Y] zeros(n+1)];
%
% Compute the model coefficients.
%
    [U,S,V]        = svd(W);
    Sdiag          = diag(S);
    indeces        = Sdiag < tol_svd;
    Sdiag(indeces) = tol_svd;
    Sinv           = diag(1./Sdiag);
    lambda         = V * Sinv *U'* b;
%
% Retrieve the model coefficients.
%
    g = lambda(m+2:m+n+1);
    for j = 1:m
       H = H + lambda(j)*(Y(:,j)*Y(:,j)');
    end
else
%
% Compute a complete quadratic model by solving a minimum least squares problem.
%
    phi_Q = zeros(m,n*(n+1)/2);
    mask  = tril(true(n));
    for i = 1:m
        y           = Y(:,i);
        aux_H       = y*y'-1/2*diag(y.^2);
        aux         = aux_H(mask);
        phi_Q(i,:)  = aux';
    end
    W = [ones(m,1) Y' phi_Q];
    b = Y_values';
%
% Compute the coefficients of the model.
%
    [U,S,V]        = svd(W,0);
    Sdiag          = diag(S);
    indeces        = Sdiag < tol_svd;
    Sdiag(indeces) = tol_svd;
    Sinv           = diag(1./Sdiag);
    lambda         = V * Sinv' *U'* b;
%
% Retrieve the model coefficients.
%
    g    = lambda(2:n+1);
    cont = n+1;
    for j = 1:n
       H(j:n,j) = lambda(cont+1:cont+n-(j-1));
       cont     = cont + n - (j-1);
    end
    H = H + H' - diag(diag(H));
end
%
% End of quad_Frob.