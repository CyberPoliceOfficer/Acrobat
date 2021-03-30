clear;
tic
%% Description
% This code generates horizontal slices of the workspace on a arbitrary z,
% and computes the performance metrics. It assumes a constant orientation.

%% Physical parameters of the platform in meters and radians.

r_b = 0.052566; %Radious of the base
r_p = 0.048139; %Radious of the platform

d_b = 0.01559; %Lenght between base's anchors
d_p = 0.00800; %Lenght between platform's anchors

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

beta_k = n*(2/3)*pi + (-1).^(k)*pi;
phi_k = (-1).^(k+1)*phi;

    
%% Discretize the 3D space

workspace_max = d+h-sin(pi/3)*r_b;

dim_y = 1000;
dim_z = dim_y;

y_min = -workspace_max;
y_max = workspace_max;
dy = (y_max - y_min)/(dim_y-1);

z_min = h;
z_max = d+h;
dz = (z_max - z_min)/(dim_z-1);

%% Prellocate the binary 3D matrix

% Manipulability Index
MI = zeros(dim_y,dim_z);

% Conditioning Index
CI = zeros(dim_y,dim_z);
CIT = zeros(dim_y,dim_z);
CIR = zeros(dim_y,dim_z);

% Payload Index
PI = zeros(dim_y,dim_z);

% Accuracy Index
AI_T = zeros(dim_y,dim_z);
AI_R = zeros(dim_y,dim_z);

% Area
Area = 0;

%% Find z

x = 0;
R = eye(3);

%% Generate the workspace
for iy = 1:dim_y
    for iz = 1:dim_z
        T = [x; y_min + dy*iy; z_min + dz*iz];
        workspace = 1;
        
        alpha_k = SixRSS_inverse_kinematics (T, R, p_k, b_k, beta_k, phi_k, d, h);

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
            Area = Area + dy*dz;

            %% Compute the Manipulability Index
            MI(iy,iz) = sqrt(det(J'*J));

            %% Compute the Conditioning Index
            CI(iy,iz) = S(end)/S(1);
            St = svd(Jt);
            CIT(iy,iz) = St(end)/St(1);
            Sr = svd(Jr);
            CIR(iy,iz) = Sr(end)/Sr(1);

            %% Compute the Payload Index
            PI(iy,iz) = S(end);
            
            %% Compute the Accuracy Index
            AI_T(iy,iz) = norm(Jt,Inf);
            AI_R(iy,iz) = norm(Jr,Inf);

        end

    end
end

toc

%% Plot everything
font=12;
y = y_min:dy:y_max;
z = z_min:dz:z_max;
[Y,Z] = meshgrid(y,z);

figure
s=pcolor(Y,Z,CI');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Dexterity')
xlabel('y (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(Y,Z,CIT');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Dexterity (T)')
xlabel('y (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(Y,Z,CIR');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Dexterity (R)')
xlabel('y (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(Y,Z,AI_T');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Displacement Sensitivity (T)')
xlabel('y (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);


figure
s=pcolor(Y,Z,AI_R');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Displacement Sensitivity (R)')
xlabel('y (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);

figure
s=pcolor(Y,Z,MI');
s.LineStyle = 'none';
pbaspect([1 1 1])
title('Manipulability')
xlabel('y (meter)') 
ylabel('z (meter)') 
colorbar
set(gca,'FontSize', font);