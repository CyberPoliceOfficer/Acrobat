function J = GetJacobian (T, R, b_k, H_k, m_k, u_k)
    M_k = R*m_k + T;
    Jx = [(M_k-H_k)' cross(M_k-T, M_k-H_k)'];
    Jq = eye(6).*dot(cross(H_k-b_k,M_k - H_k),u_k);
    J = inv(Jx)*Jq;
end
