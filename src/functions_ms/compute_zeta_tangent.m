function zeta_tangent = compute_zeta_tangent(gammdot0, V0, a, sigma_n, p, tauf, f0, b, theta, Dc_prime, eta, yy_gammadot)
    
    s = eta*trapz(yy_gammadot, dgammadotdtauf(gammdot0, V0, a, sigma_n-p, tauf, f0, b, theta, Dc_prime), 2);
    
    % Compute tagent value
    zeta_tangent = 1 + s;
end