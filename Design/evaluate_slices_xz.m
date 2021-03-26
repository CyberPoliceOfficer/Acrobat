clear;
tic
%% Description
% This code generates horizontal slices of the workspace on a arbitrary z,
% and computes the performance metrics. It assumes a constant orientation.

%% Physical parameters of the platform in meters and radians.

r_b = 0.052566; %Radious of the base
r_p = 0.048139; %Radious of the platform

d_b = 0.01559; %Lenght between base's anchors
d_p = 0.0800; %Lenght between platform's anchors

h = 0.027; %Servo's arm lenght
d = 0.1175; %Platform's arm lenght

phi = 20*2*pi/360; %Angle between servo's arm and platform's base

%% Compute vector bk and pk

k = [1 2 3 4 5 6];
n = floor((k-1)/2);

theta_b = n*(2/3)*pi + (-1).^k*asin(d_b/(2*r_b));
theta_p = n*(2/3)*pi + (-1).^k*asin(d_p/(2*r_p));

b_k = [r_b * cos(theta_b); r_b * sin(theta_b);zeros(1,6)];
p_k = [r_p * cos(theta_p); r_p * sin(theta_p);zeros(1,6)];


%% Compute beta and gamma

beta_k = n*(2/3)*pi + (-1).^(k)*pi/2;
phi_k = (-1).^(k+1)*phi;

    
%% Discretize the 3D space

workspace_max = d+h-sin(pi/3)*r_b;

dim_x = 1000;
dim_z = dim_x;

x_min = -workspace_max;
x_max = workspace_max;
dx = (x_max - x_min)/(dim_x-1);

z_min = h;
z_max = d+h;
dz = (z_max - z_min)/(dim_z-1);

%% Prellocate the binary 3D matrix

% Manipulability Index
MI = zeros(dim_x,dim_z);

% Conditioning Index
CI = zeros(dim_x,dim_z);
CIT = zeros(dim_x,dim_z);
CIR = zeros(dim_x,dim_z);

% Payload Index
PI = zeros(dim_x,dim_z);

% Accuracy Index
AI_T = zeros(dim_x,dim_z);
AI_R = zeros(dim_x,dim_z);

% Area
Area = 0;

%% Find z

y = 0;

%% Generate the workspace
for ix = 1:dim_x
    for iz = 1:dim_z
        T = [x_min + dx*ix; y; z_min + dz*iz];
        workspace = 1;
        
        alpha_k = SixRSS_inverse_kinematics (T, eye(3), p_k, b_k, beta_k, phi_k, d, h);

        if (isreal(alpha_k)== 0)
           workspace = 0;
        end

        if (workspace)
            %% Pre compute a few things
            h_k = h* [sin(beta_k).*sin(phi_k).*sin(alpha_k) + cos(beta_k).*cos(alpha_k);
                      -cos(beta_k).*sin(phi_k).*sin(alpha_k) + sin(beta_k).*cos(alpha_k);
                      cos(phi_k).*sin(alpha_k)];

            H_k = b_k + h_k;
            u_k = -[cos(beta_k+pi/2).*cos(phi_k);
                   sin(beta_k+pi/2).*cos(phi_k);
                   sin(phi_k)];
            J = GetJacobian (T, eye(3), b_k, H_k, p_k, u_k);
            [~,S,~] = svd(J);
            S = diag(S);
            Det_J = abs(prod(S));
            Jr = J(4:6,:);
            Jt = J(1:3,:);

            %% Compute the Volume
            Area = Area + dx*dz;

            %% Compute the Manipulability Index
            MI(ix,iz) = sqrt(det(J'*J));

            %% Compute the Conditioning Index
            CI(ix,iz) = S(end)/S(1);
            St = svd(Jt);
            CIT(ix,iz) = St(end)/St(1);
            Sr = svd(Jr);
            CIR(ix,iz) = Sr(end)/Sr(1);

            %% Compute the Payload Index
            PI(ix,iz) = S(end);
            
            %% Compute the Accuracy Index
            AI_T(ix,iz) = norm(Jt,Inf);
            AI_R(ix,iz) = norm(Jr,Inf);

        end

    end
end

toc

%% Plot everything
font=12;
x = x_min:dx:x_max;
z = z_min:dz:z_max;
[X,Z] = meshgrid(x,z);

figure
s=pcolor(X,Z,CI');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Dexterity')
xlabel('x (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(X,Z,CIT');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Dexterity (T)')
xlabel('x (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(X,Z,CIR');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Dexterity (R)')
xlabel('x (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(X,Z,AI_T');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Displacement Sensitivity (T)')
xlabel('x (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);


figure
s=pcolor(X,Z,AI_R');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Displacement Sensitivity (R)')
xlabel('x (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(X,Z,MI');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Manipulability')
xlabel('x (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);