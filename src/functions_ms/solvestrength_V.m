function [tauf_return, nb_it_return] = solvestrength_V(tauf, sigma_n, p, theta, tau_dyn, V, gammadot0, a, f0, b, Dc_prime, eta, yy_gammadot)

    % Newton-Raphson parameters
    tol = 1e-7;          % Convergence tolerance
    maxit = 200;          % Maximum number of iterations
    nb_it = 0;           % Iteration counter
    d_tauf = Inf;        % Initialize change in tauf

    % Newton-Raphson loop to solve for shear strength tauf
    while (abs(d_tauf) > tol) && (nb_it < maxit)
        % Residual: difference between applied and constitutive shear stress
        d_zeta = -compute_zeta_V(tauf, eta, V) + tau_dyn;
       
        % Derivative of residual w.r.t. tauf
        zeta_p = compute_zeta_tangent_V(gammadot0, a, sigma_n, p, tauf, f0, b, theta, Dc_prime, eta, yy_gammadot);
        
        % Newton-Raphson update
        d_tauf = d_zeta / zeta_p;
        tauf = tauf + d_tauf;
        nb_it = nb_it + 1;
    end

    % Check convergence
    if nb_it == maxit
        fprintf("Newton-Raphson failed at finding an equilibrated solution (error %.5f)\n", abs(d_tauf));
    end

    % Return the updated shear strength and iteration count
    tauf_return = tauf;
    nb_it_return = nb_it;
end
