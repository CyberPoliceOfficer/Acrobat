function rot = GetRotMatv0 (gamma, theta, psi)
    rot = Rx(gamma)*Ry(theta)*Rz(psi);  
end
