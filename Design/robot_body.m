clear;
K1 = 1e-03;
K2 = 6.4e-03*K1;
a = sqrt(2);

factor = 0;

%% Unnamed Robot
theta = deg2rad([0 180 90 90 90 90]);
gamma = deg2rad([0 0 0 180 90 270]);
alpha = deg2rad([90 90 0 180 90 90]);
beta = deg2rad([90 270 90 270 180 0]);
r = [cos(gamma).*sin(theta); sin(gamma).*sin(theta); cos(theta)];
u = [cos(beta).*sin(alpha); sin(beta).*sin(alpha); cos(alpha)];

w = [-1 1 -1 1 -1 1];
A = [u; cross(r,u) - factor*w.*u];
S = svd(A);
kinv_Ubot=S(6)/S(1);

Ainv = inv(A);
b_Ubot = Ainv(:,1:3);
c_Ubot = Ainv(:,4:6);


figure
hold on
axis equal
O = zeros(3,6);
for i=1:6
    text(r(1,i),r(2,i),r(3,i),string(i),'FontSize',14)
    v=[O(:,i)';r(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'blue')
end

ul = r + u;

for i=1:6
    v=[r(:,i)';ul(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'red')
end

%% Acrobat
theta = deg2rad([90 210 330 60 180 300]);
phi = deg2rad([0 120 240 60 180 300]);
gamma = deg2rad([90 90 90 45 45 45]);
z = [1 1 1 -1 -1 -1];
w = [1 1 1 1 1 1];


theta = deg2rad([0 330 120 90 240 210]);
phi = deg2rad([0 60 120 180 240 300]);
gamma = deg2rad([35.26439 90 35.26439 90 35.26439 90]);

z = [-1 1 -1 1 -1 1];

% Plot some stuff to debug
r = [a/sqrt(3)*cos(phi); a/sqrt(3)*sin(phi); z*a/sqrt(6)];
u = [cos(theta).*sin(gamma); sin(theta).*sin(gamma); cos(gamma)];

% Compute matrix manually
A_acrobat = [K1*cos(theta).*sin(gamma); 
    K1*sin(theta).*sin(gamma);
    K1*cos(gamma);
    K1*(a/sqrt(3)*sin(phi).*cos(gamma) - a/sqrt(6)*z.*sin(theta).*sin(gamma)) - K2*w.*cos(theta).*sin(gamma);
    K1*(a/sqrt(6)*z.*cos(theta).*sin(gamma) - a/sqrt(3)*cos(phi).*cos(gamma)) - K2*w.*sin(theta).*sin(gamma);
    K1*(a/sqrt(3)*cos(phi).*sin(theta).*sin(gamma) - a/sqrt(3)*sin(phi).*cos(theta).*sin(gamma)) - K2*w.*cos(gamma)];

% Use matlab
A_acrobat = [u; cross(r,u) - factor*w.*u];
S = svd(A_acrobat);
kinv_acrobat=S(6)/S(1);
            
            
A = inv(A_acrobat);
b_acrobat = A(:,1:3);
c_acrobat = A(:,4:6);



figure
hold on
axis equal
O = zeros(3,6);
for i=1:6
    text(r(1,i),r(2,i),r(3,i),string(i),'FontSize',14)
    v=[O(:,i)';r(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'blue')
end

ul = r + u;

for i=1:6
    v=[r(:,i)';ul(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'red')
end

%% Space CoBot
theta = deg2rad([0 60 120 180 240 300]);
phi = deg2rad([55 -55 55 -55 55 -55]);
w = [-1 1 -1 1 -1 1];
d=1;

A_spacecobot = [1*sin(theta).*sin(phi); 
    -1*cos(theta).*sin(phi);
    1*cos(phi);
    (1*d*cos(phi) - factor*w.*sin(phi)).*sin(theta);
    -(1*d*cos(phi) - factor*w.*sin(phi)).*cos(theta);
    -1*d*sin(phi) - factor*w.*cos(phi)];
S = svd(A_spacecobot);
kinv_spacecobot=S(6)/S(1);


A = inv(A_spacecobot);
b_spacecobot = A(:,1:3);
c_spacecobot = A(:,4:6);

%% Validation

F_max_acrobat = min(1./vecnorm(b_acrobat'));
M_max_acrobat = min(1./vecnorm(c_acrobat'));

F_max_spacecobot = min(1./vecnorm(b_spacecobot'));
M_max_spacecobot = min(1./vecnorm(c_spacecobot'));

F_max_Ubot = min(1./vecnorm(b_Ubot'));
M_max_Ubot = min(1./vecnorm(c_Ubot'));