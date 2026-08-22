function zeta_tangent = compute_zeta_tangent_V(gammadot0, a, sigma_n, p, tauf, f0, b, theta, Dc_prime, eta, yy_gammadot)
    
    s = eta*trapz(yy_gammadot, dVdtauf(gammadot0, a, sigma_n-p, tauf, f0, b, theta, Dc_prime), 2);
    
    % Compute tagent value
    zeta_tangent = 1 + s;
end