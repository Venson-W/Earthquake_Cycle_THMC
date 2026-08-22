function zeta = compute_zeta(V0, a, sigma_n, p, tauf, f0, b, theta, Dc, eta)

    s = eta .* shear_rate(tauf, theta, V0, f0, b, a, Dc, sigma_n, p);

    zeta = tauf + s;
end