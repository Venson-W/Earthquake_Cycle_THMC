function zeta = compute_zeta(gammadot0, V0, a, sigma_n, p, tauf, f0, b, theta, Dc_prime, eta, yy_gammadot)

    s = eta .* trapz(yy_gammadot, shear_strain(tauf, theta, gammadot0, V0, f0, b, a, Dc_prime, sigma_n, p), 2);

    zeta = tauf + s;
end
