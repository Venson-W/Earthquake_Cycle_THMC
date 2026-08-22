function zeta_tangent = compute_zeta_tangent(V0, a, sigma_n, p, tauf, f0, b, theta, Dc, eta)
    
    s = eta* dVdtauf(V0, a, sigma_n-p, tauf, f0, b, theta, Dc);
    
    % Compute tagent value
    zeta_tangent = 1 + s;
end