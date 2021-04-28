function alpha_k = SixRSS_inverse_kinematics (T, R, p_k, b_k, beta_k, phi_k, d, h)
    %% Compute vector ik

    i_k = repmat (T, [1,6]) + R*p_k - b_k;
    

    %% Compute alpha

    f_k = (cos(beta_k).*i_k(1,:) + sin(beta_k).*i_k(2,:))*2*h;

    e_k = (+sin(beta_k).*sin(phi_k).*i_k(1,:) - cos(beta_k).*sin(phi_k).*i_k(2,:) + cos(phi_k).*i_k(3,:))*2*h;

    g_k = vecnorm(i_k).^2 -(d^2-h^2);

    alpha_k = asin(g_k./sqrt(e_k.^2 + f_k.^2)) - atan2(f_k,e_k);
end