clear;

x = [0.05257; 0.04814; 0.1624; 0.0832; 0.1175; 0.027;0.3491; pi/2];

%% Desired translation and rotation of the platform

%T = [0 -0.0808578261898999 0]';%Translation
%T = [0; 0; 0.11];
%T = [0; 0; 0.1];
W = [0; 0; 0]; %Rotation
%R = GetRotMatv0(W(1),W(2),W(3));

%% Physical parameters of the platform in meters and radians.
r_b = x(1); r_p = x(2); d_b = x(3); d_p = x(4); d = x(5); h = x(6); phi = x(7); beta = x(8);

%% Compute vector bk and pk

k = [1 2 3 4 5 6];
n = floor((k-1)/2);

theta_b = n*(2/3)*pi + (-1).^k*d_b;
theta_p = n*(2/3)*pi + (-1).^k*d_p;

b_k = [r_b * cos(theta_b); r_b * sin(theta_b);zeros(1,6)];
p_k = [r_p * cos(theta_p); r_p * sin(theta_p);zeros(1,6)];


%% Compute beta and gamma

beta_k = n*(2/3)*pi + (-1).^(k)*beta;
phi_k = (-1).^(k+1)*phi;

%% Home
home = mean((d^2-(r_p*cos(theta_p) - r_b*cos(theta_b) - h*cos(beta_k)).^2 - (r_p*sin(theta_p) - r_b*sin(theta_b) - h*sin(beta_k)).^2).^(1/2))

T = [0; 0; home];
R = eye(3);

%% Compute the inverse kinematics

alpha_k = SixRSS_inverse_kinematics (T, R, p_k, b_k, beta_k, phi_k, d, h);

%% Check if the restriction still stands

h_k = h* [sin(beta_k).*sin(phi_k).*sin(alpha_k) + cos(beta_k).*cos(alpha_k);
          -cos(beta_k).*sin(phi_k).*sin(alpha_k) + sin(beta_k).*cos(alpha_k);
          cos(phi_k).*sin(alpha_k)];
H_k = b_k + h_k;

u_k = -[cos(beta_k+pi/2).*cos(phi_k);
       sin(beta_k+pi/2).*cos(phi_k);
       sin(phi_k)];
           
P_k = repmat (T, [1,6]) + R*p_k;
norm_d_k = vecnorm(P_k - H_k);
norm_h_k = vecnorm(h_k);

alpha_k

J = GetJacobian (T, R, b_k, H_k, p_k, u_k);
M_k = R*p_k + T;
Jx = [(M_k-H_k)' cross(M_k-T, M_k-H_k)'];
Jq = eye(6).*dot(cross(H_k-b_k,M_k - H_k),u_k);
 
[~,S,~] = svd(J);
S = diag(S);
% prod(S)
Det_J = det(J);
GCI= S(end)/S(1);
Jr = J(4:6,:);
Jt = J(1:3,:);
GCI_2 = 1/(norm(J,2)* norm(inv(J),2));
sigma_t = norm(Jt,Inf);
sigma_r = norm(Jr,Inf);


%% Draw everything using alpha
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
    v=[b_k(:,i)';b_k(:,i)'+u_k(:,i)'];
    text(b_k(1,i),b_k(2,i),b_k(3,i),string(i),'FontSize',14)
    %plot3(v(:,1),v(:,2),v(:,3),'green')
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