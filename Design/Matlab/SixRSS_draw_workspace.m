clear;
tic

%% Physical parameters of the platform in meters and radians.

r_b = 0.5; %Radious of the base
r_p = 0.3; %Radious of the platform

d_b = 0.2; %Lenght between base's anchors
d_p = 0.4; %Lenght between platform's anchors

h = 0.3; %Servo's arm lenght
d = 0.7; %Platform's arm lenght

phi = 0.3491; %Angle between servo's arm and platform's base



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

%% Compute the binary 3D matrix

%% Discretize the 3D space

Nx=50;
Ny=50;
Nz=50;
workspace_max = h+d - cos(pi/3)*r_b;
x = linspace(-workspace_max,workspace_max,Nx);
y = linspace(-workspace_max,workspace_max,Ny);
z = linspace(0,d+h,Nz);
[X,Y,Z] = meshgrid(x,y,z);
density = (2*workspace_max)/(Nx-1)*(2*workspace_max)/(Ny-1)*(d-h)/(Nz-1);

dim_x = size(X); dim_x = dim_x(1);
dim_y = size(Y); dim_y = dim_y(2);
dim_z = size(Z); dim_z = dim_z(3);

workspace = zeros(dim_x,dim_y,dim_z);

%% Generate the workspace
for ix = 1:dim_x
    for iy = 1:dim_y
        for iz = 1:dim_z
            T = [X(ix,iy,iz); Y(ix,iy,iz); Z(ix,iy,iz)];
            alpha_k = SixRSS_inverse_kinematics (T, eye(3), p_k, b_k, beta_k, phi_k, d, h);

            if (isreal(alpha_k)== 1)
               workspace(ix,iy,iz) = 1;
            end
              
        end
    end
end

toc

volume = sum(workspace(:))*density;

if 0
%% Draw the workspace
figure
p = patch(isosurface(X,Y,Z,workspace,0));
alpha(0.5);
hold on;
isonormals(X,Y,Z,workspace,p)
daspect([1 1 1])
view(3);
axis tight
camlight('headlight', 'infinite') 
lighting gouraud

%% Desired translation and rotation of the platform
z = sqrt(d^2 - (r_p*cos(theta_p) - r_b*cos(theta_b) + h*cos(beta_k)).^2 - (r_p*sin(theta_p) - r_b*sin(theta_b) + h*sin(beta_k)).^2);
   
T = [0; 0; mean(z)]; %Translation

R = eye(3);

%% Compute the inverse kinematics

alpha_k = SixRSS_inverse_kinematics (T, R, p_k, b_k, beta_k, phi_k, d, h);

%% Check if the restriction still stands

h_k = h* [sin(beta_k).*sin(phi_k).*sin(alpha_k) + cos(beta_k).*cos(alpha_k);
                      -cos(beta_k).*sin(phi_k).*sin(alpha_k) + sin(beta_k).*cos(alpha_k);
                      cos(phi_k).*sin(alpha_k)];

H_k = b_k + h_k;
           
P_k = repmat (T, [1,6]) + R*p_k;


%% Draw everything using alpha
patch(b_k(1,:), b_k(2,:), b_k(3,:),'green');
daspect([1 1 1])
patch(P_k(1,:), P_k(2,:), P_k(3,:),'red');
alpha(0.5);
set(gca,'FontSize', 18);


for i=1:6
    v=[H_k(:,i)';b_k(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'blue')
end

for i=1:6
    v=[H_k(:,i)';P_k(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'blue')
end
 
%OptionZ.FrameRate=30;OptionZ.Duration=6;OptionZ.Periodic=true;
%CaptureFigVid([-20,20;-110,20;-190,20;-290,20;-380,20], 'WellMadeVid',OptionZ)
end
