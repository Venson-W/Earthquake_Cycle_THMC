function res = residual_tauf(tauf, k, Vp, t, x, tau0, eta, gammdot0, a, b, f0, Dc_prime, thetan, dt, sigma_n, p, yy_gammadot)
    % Effective normal stress
    sigma_eff = sigma_n - p;

    % --- Step 1: Estimate slip velocity V1 using thetan ---
    gammadot1 = 2 .* gammdot0 .* sinh(tauf ./ (a .* sigma_eff)) .* ...
         exp(-1./a.*(f0+b*log(thetan*gammdot0/Dc_prime)));
   
    % --- Step 2: Update state variable theta using aging law ---
    theta = state_theta(thetan, gammadot1, dt, Dc_prime);

    % --- Step 3: Recompute slip velocity V using updated theta ---
    gammadot = 2 .* gammdot0 .* sinh(tauf ./ (a .* sigma_eff)) .* ...
         exp(-1./a.*(f0+b*log(theta*gammdot0/Dc_prime)));

    % --- Step 4: Compute residual ---
    res = gather(k * (Vp * t - x) - tauf + tau0 - eta .* slip_velocity(gammadot,yy_gammadot));
end
