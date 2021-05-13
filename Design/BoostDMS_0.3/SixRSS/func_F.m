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
% 1 = r_b, 2 = r_m, 3 = d_b, 4 = d_m, 5 = d, 6 = h, 7 = phi, 8 = beta
%F = EvaluateWorkspaceDiscretization(x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8));
r_b = x(1); r_m = x(2); d_b = x(3); d_m = x(4); d = x(5); h = x(6); phi = x(7); beta = x(8);

%% Options
Nx = 50;
Ny = 50;
Nz = 50;
n_metrics = 3;
joint_limit = 1;

%% Compute vector bk and pk
k = 1:6;
n = floor((k-1)/2);

theta_b = n*(2/3)*pi + (-1).^k*d_b;
theta_p = n*(2/3)*pi + (-1).^k*d_m;

b_k = [r_b * cos(theta_b); r_b * sin(theta_b);zeros(1,6)];
m_k = [r_m * cos(theta_p); r_m * sin(theta_p);zeros(1,6)];

%% Compute beta and gamma
beta_k = n*(2/3)*pi + (-1).^(k)*beta;
phi_k = (-1).^(k+1)*phi;

%% Compute u_k
u_k = -[cos(beta_k+pi/2).*cos(phi_k);
        sin(beta_k+pi/2).*cos(phi_k);
        sin(phi_k)];

%% Discretize the 3D space
workspace_max_xy = h+d - cos(pi/3)*r_b;
workspace_min_z = 0;
workspace_min_xy = -workspace_max_xy; 
dx = (2*workspace_max_xy)/(Nx-1);
dy = (2*workspace_max_xy)/(Ny-1);
dz = (d+h)/(Nz-1);
density = dx*dy*dz;

F = zeros(n_metrics,1);
n_nodes = 0;
GTSI = 0;
GRSI = 0;

R = eye(3);

%% Evaluate the workspace
for ix = 1:Nx
    for iy = 1:Ny
        for iz = 1:Nz
            T = [workspace_min_xy + (ix-1)*dx; workspace_min_xy + (iy-1)*dy; workspace_min_z + (iz-1)*dz];
            i_k = repmat (T, [1,6]) + R*m_k - b_k;

            %% Kinematic constraints

            f_k = (cos(beta_k).*i_k(1,:) + sin(beta_k).*i_k(2,:))*2*h;

            e_k = (+sin(beta_k).*sin(phi_k).*i_k(1,:) - cos(beta_k).*sin(phi_k).*i_k(2,:) + cos(phi_k).*i_k(3,:))*2*h;

            g_k = vecnorm(i_k).^2 -(d^2-h^2);

            t_k = sqrt(e_k.^2 + f_k.^2);

            % Check if inside the workspace
            if (any(abs(g_k) > t_k))
                continue;
            end

            alpha_k = asin(g_k/t_k) - atan2(f_k, e_k);
            h_k = h* [sin(beta_k).*sin(phi_k).*sin(alpha_k) + cos(beta_k).*cos(alpha_k);
                      -cos(beta_k).*sin(phi_k).*sin(alpha_k) + sin(beta_k).*cos(alpha_k);
                      cos(phi_k).*sin(alpha_k)];
            MR_k = R*m_k;
            d_k = MR_k + T - h_k - b_k;

            % Check if first joint constraint is violated
            %j1_k = cross(h_k, d_k);
            %if(any(atan2(vecnorm(cross(j1_k,u_k)),dot(j1_k,u_k)) > joint_limit))
            %    continue;
            %end
            %if(any(abs(atan2(vecnorm(cross(d_k,MR_k)),dot(d_k,MR_k)) - pi/2) > joint_limit))
            %    continue;
            %end
            
            Jq = dot(cross(h_k, d_k),u_k);
            InvJ = [(d_k./Jq)' (cross(MR_k, d_k)./Jq)'];
            J = inv(InvJ);

            % Compute the metrics
            n_nodes = n_nodes + 1;
            GTSI = GTSI + norm(J(1:3,:),Inf);
            GRSI = GRSI + norm(J(4:6,:),Inf);

        end
    end
end

if (n_nodes > 0)
    F(2) = GTSI/n_nodes;
    F(3) = GRSI/n_nodes;
    F(1) = -n_nodes*density;
end
%-----------------------------------------------------------------------------
%  End of block to be user modified.
% -----------------------------------------------------------------------------
%
% End of func_F.