function [xi_next, Jn] = reaction_progress(xi_n, Tn, T, dt, A, Q, R)
%STEP_XI_ARRHENIUS_TN  One-step update for xi'=(1-xi)*A*exp(-Q/(R*T)).
%   T is a spatial field known only at the previous time step (frozen over dt).
%
%   Inputs:
%     xi_n : current xi field (vector/matrix/ND array)
%     Tn   : temperature field at t_n, same size as xi_n (>0 everywhere)
%     dt   : time step (scalar, >0)
%     A,Q,R: Arrhenius parameters (scalars or arrays broadcastable to xi_n)
%
%   Outputs:
%     xi_next : xi at t_{n+1}
%     Jn      : per-location step integral ∫ k dt = dt*A*exp(-Q/(R*Tn))
%
%   Notes:
%     • Exact for the reaction with frozen T over the step.
%     • Unconditionally stable, preserves 0 ≤ xi ≤ 1 (up to roundoff).

    % basic checks
    if ~isscalar(dt) || dt <= 0
        error('dt must be a positive scalar.');
    end
    if ~isequal(size(Tn), size(xi_n))
        error('Tn must have the same size as xi_n.');
    end

    % Arrhenius rate at t_n (elementwise, supports arrays)
    a  = Q ./ R;                        % works for scalar or array Q,R
    kN = A .* exp( - a ./ Tn );         % k(t_n, x)
    kNpone = A .* exp( - a ./ T );
    % Step integral and exponential update (use expm1 for small Jn)
    Jn      = 0.5 * dt .* (kN + kNpone);                 % ∫_{t_n}^{t_{n+1}} k ≈ dt*k(t_n)
    xi_next = xi_n - (1 - xi_n) .* expm1(-Jn);

    % clip tiny numerical drift
    xi_next = max(0, min(1, xi_next));
end
