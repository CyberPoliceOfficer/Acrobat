clear;
%% Desired translation and rotation of the platform

T = [0.00 0.02 0.6]';%Translation
W = [0; 0; 0]; %Rotation
R = GetRotMatv0(W(1),W(2),W(3));

%% Physical parameters of the platform in meters and radians.

r_b = 0.5; %Radious of the base
r_p = 0.3; %Radious of the platform

d_b = 0.2; %Lenght between base's anchors
d_p = 0.4; %Lenght between platform's anchors

d = 0.7; %Platform's arm lenght
h = 0.3; %Servo's arm lenght

phi = 0; %Angle between servo's arm and platform's base


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

%% Compute the inverse kinematics

alpha_k = SixRSS_inverse_kinematics (T, R, p_k, b_k, beta_k, phi_k, d, h);


%% Check if the restriction still stands

home = (d^2-(r_p*cos(theta_p) - r_b*cos(theta_b) + h*cos(beta_k)).^2 - (r_p*sin(theta_p) - r_b*sin(theta_b) + h*sin(beta_k)).^2).^(1/2);

h_k = h* [sin(beta_k).*sin(phi_k).*sin(alpha_k) + cos(beta_k).*cos(alpha_k);
          -cos(beta_k).*sin(phi_k).*sin(alpha_k) + sin(beta_k).*cos(alpha_k);
          cos(phi_k).*sin(alpha_k)];
H_k = b_k + h_k;

u_k = -[cos(beta_k+pi/2).*cos(phi_k);
       sin(beta_k+pi/2).*cos(phi_k);
       sin(phi_k)];
   

%% Draw everything using alpha
P_k = repmat (T, [1,6]) + R*p_k;
figure;
patch(b_k(1,:), b_k(2,:), b_k(3,:),'black');
hold on;
daspect([1 1 1])
patch(P_k(1,:), P_k(2,:), P_k(3,:),'red');

I_k = repmat (T, [1,6]) + R*p_k;

if false
    for i=1:6
        v=[I_k(:,i)';b_k(:,i)'];
        plot3(v(:,1),v(:,2),v(:,3),'red')
    end
end

for i=1:6
    v=[H_k(:,i)';b_k(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'blue')
end

for i=1:6
    v=[H_k(:,i)';P_k(:,i)'];
    plot3(v(:,1),v(:,2),v(:,3),'blue')
end

if false
    for i=1:6
        v=[b_k(:,i)';u_k(:,i)'];
        plot3(v(:,1),v(:,2),v(:,3),'green')
    end
end

%% Start conditions
x_k = [0 0 mean(home) 0 0 0]';
n=1;
%% Inverse kinematics using Newton Rhapson method
while (1)
    
    x_kminus = x_k;
    T = x_k(1:3);
    R = GetRotMatv0(x_k(4),x_k(5),x_k(6));
    
    alpha = SixRSS_inverse_kinematics (T, R, p_k, b_k, beta_k, phi_k, d, h);
    h_k = h* [sin(beta_k).*sin(phi_k).*sin(alpha) + cos(beta_k).*cos(alpha);
          -cos(beta_k).*sin(phi_k).*sin(alpha) + sin(beta_k).*cos(alpha);
          cos(phi_k).*sin(alpha)];
    H_k = b_k + h_k;
    J = GetJacobian (T, R, b_k, H_k, p_k, u_k);
      
    P_k = repmat (T, [1,6]) + R*p_k;
    f_i = vecnorm(P_k - H_k,2)- d;
    
    %x_k = x_kminus + mldivide(J,-f_i');
    x_k = x_kminus + J*(alpha_k - alpha)';
    n=n+1;
    if(norm(x_k - x_kminus)<1e-4 || n >= 1e2)
       break 
    end
end

