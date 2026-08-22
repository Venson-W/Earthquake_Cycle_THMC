function [T, p] = implicit_p_T(Tn, pn, gamma, tauf, dt, ...
        ~, ~, ...       % T_boundarynorth, T_boundarysouth  (unused: Robin BC)
        p_undrained_poroelas_N, p_undrained_poroelas_S, ...       % p_boundarynorth, p_boundarysouth  (unused: Robin BC)
        M_alpha_hy, alpha_th, rho_c, Lambda, dy, ...
        n_pl_term, use_plastic_term, decomposition, ...
        rect_m_frac, Pr, Er, alpha_hy_bulk, alpha_th_bulk)
%IMPLICIT_P_T  Conservative staggered-FD solver for fault-normal p and T
%              with diffusive (Robin) boundary conditions at both gouge edges.
%
% STAGGERED GRID (n physical nodes, n+1 original faces):
%
%   [phys_1] -f2- [phys_2] - ... - [phys_{n-1}] -f_n- [phys_n]
%      |                                                    |
%    f_1 (south boundary face)              f_{n+1} (north boundary face)
%
% GHOST-NODE EXTENSION (n+2 nodes, n+3 faces):
%
%   [ghost_S] -f1*- [phys_1] -f2- ... -f_n- [phys_n] -f_{n+2}*- [ghost_N]
%      |                                                               |
%   f0* (outer south, Robin)                          f_{n+3}* (outer north, Robin)
%
% Face index mapping (1-based, n+3 total):
%   k = 1      : outer south Robin face,  conductance = sigma * dy
%   k = 2..n+2 : gouge faces shifted from M_alpha_hy(1..n+1)
%   k = n+3    : outer north Robin face,  conductance = sigma * dy
%
% Ghost nodes are quasi-stationary (no d/dt). Their rows encode the Robin
% boundary condition:
%
%   M_face * (p_phys - p_ghost) / dy = M_outer * p_ghost / dy
%   => p_ghost = M_face / (M_outer + M_face) * p_phys
%
% with M_outer = sigma * dy = sqrt(alpha_bulk / (pi*dt)) * dy.
%
% When sigma -> inf  the ghost -> 0  (Dirichlet, fast bulk diffusion).
% When sigma -> 0    the ghost -> p_phys  (Neumann, impermeable bulk).
%
% T_boundarynorth/south and p_boundarynorth/south are retained for API
% compatibility with the calling code but are not used.

    n = length(Tn);    % number of physical nodes
    N = n + 2;         % extended system size (physical + 2 ghost nodes)

    if n < 2
        error('implicit_p_T: need at least 2 physical nodes.');
    end

    Tn = Tn(:);
    pn = pn(:);

    % ── Robin half-space impedances  σ = sqrt(α_bulk / (π·dt))  [m/s] ──
    sigma_T = sqrt(alpha_th_bulk / (pi * dt));
    sigma_p = sqrt(alpha_hy_bulk / (pi * dt));

    % ── Build extended face arrays (n+3 faces) ───────────────────────────
    %
    % M_alpha_hy has n+1 elements: M(1) = south boundary face (between
    % ghost_S and phys_1 in the extended system), M(2..n) = interior faces,
    % M(n+1) = north boundary face (between phys_n and ghost_N).
    % Shift these to positions 2..n+2 and add Robin outer faces at 1 and n+3.

    M_p_ext        = zeros(1, n+3);
    M_p_ext(1)     = sigma_p * dy;           % outer south Robin face
    M_p_ext(2:n+2) = M_alpha_hy(1:n+1);     % gouge faces (shifted by 1)
    M_p_ext(n+3)   = sigma_p * dy;           % outer north Robin face

    % Temperature: uniform alpha_th on all gouge faces
    M_T_ext        = zeros(1, n+3);
    M_T_ext(1)     = sigma_T * dy;           % outer south Robin face
    M_T_ext(2:n+2) = alpha_th;               % uniform gouge thermal faces
    M_T_ext(n+3)   = sigma_T * dy;           % outer north Robin face

    % ══════════════════════════════════════════════════════════════════════
    % TEMPERATURE  —  (N x N) tridiagonal
    % ══════════════════════════════════════════════════════════════════════
    %
    % Extended node layout (MATLAB 1-indexed):
    %   index 1     : ghost_S
    %   index 2..N-1: physical nodes (physical index k = i-1)
    %   index N     : ghost_N
    %
    % Implicit time step for physical node k:
    %   (T_k^{n+1} - Tn_k) / dt =
    %       [M_T_ext(i+1)*(T_{k+1} - T_k) - M_T_ext(i)*(T_k - T_{k-1})] / dy^2
    %       + source_k
    %
    % Ghost-node row (quasi-stationary, no d/dt):
    %   M_T_ext(2)*(T_phys1 - T_ghost_S) - M_T_ext(1)*T_ghost_S = 0
    %   => T_b*T_ghost + T_c*T_phys1 = 0

    T_a = zeros(N, 1);
    T_b = zeros(N, 1);
    T_c = zeros(N, 1);
    T_d = zeros(N, 1);

    % Ghost south (i = 1)
    T_b(1) =  M_T_ext(1) + M_T_ext(2);
    T_c(1) = -M_T_ext(2);
    % T_d(1) = 0  (already zero)

    % Physical nodes (i = 2 .. N-1, physical index k = i-1)
    for i = 2:(N-1)
        k = i - 1;                                 % physical index  1..n

        T_a(i) = -dt * M_T_ext(i)   / dy^2;       % south-face flux coeff
        T_c(i) = -dt * M_T_ext(i+1) / dy^2;       % north-face flux coeff
        T_b(i) =  1 - T_a(i) - T_c(i);

        % Shear heating + chemical decomposition source.
        % Physical boundary nodes (k=1, k=n) sit at the gouge–rock interface
        % where shear strain rate is negligible; their source is zero.
        if k >= 2 && k <= n-1
            gidx  = k - 1;           % source array offset (same as original)
            src_T = tauf * gamma(gidx) / rho_c ...
                    - rect_m_frac * Er * decomposition(gidx);
        else
            src_T = 0;
        end

        T_d(i) = Tn(k) + dt * src_T;
    end

    % Ghost north (i = N)
    T_a(N) = -M_T_ext(N);
    T_b(N) =  M_T_ext(N) + M_T_ext(N+1);
    % T_d(N) = 0  (already zero)

    T_ext = thomas_algorithm(T_a, T_b, T_c, T_d, N);
    T     = T_ext(2:N-1);     % extract n physical-node values

    % ══════════════════════════════════════════════════════════════════════
    % PRESSURE  —  (N x N) tridiagonal
    % ══════════════════════════════════════════════════════════════════════
    delta_T = (T - Tn) / dt;   % dT/dt for thermoelastic pressurisation

    p_a = zeros(N, 1);
    p_b = zeros(N, 1);
    p_c = zeros(N, 1);
    p_d = zeros(N, 1);

    
    % Ghost south (i = 1)
    p_b(1) =  M_p_ext(1) + M_p_ext(2);   
    p_c(1) = -M_p_ext(2);                
    p_d(1) =  M_p_ext(1) * p_undrained_poroelas_S;  

    % Physical nodes (i = 2 .. N-1, physical index k = i-1)
    for i = 2:(N-1)
        k = i - 1;                                 % physical index  1..n

        p_a(i) = -dt * M_p_ext(i)   / dy^2;
        p_c(i) = -dt * M_p_ext(i+1) / dy^2;
        p_b(i) =  1 - p_a(i) - p_c(i);

        % Thermoelastic pressurisation + chemical decomposition (all nodes)
        src_p = Lambda * delta_T(k) + rect_m_frac * Pr * decomposition(k);

        % Dilatancy term: interior physical nodes only (same gidx as T)
        if k >= 2 && k <= n-1
            gidx = k - 1;
            if use_plastic_term
                src_p = src_p + n_pl_term(gidx);
            end
        end

        p_d(i) = pn(k) + dt * src_p;
    end

    % Ghost north (i = N)
    p_a(N) = -M_p_ext(N);                 
    p_b(N) =  M_p_ext(N) + M_p_ext(N+1);
    p_d(N) =  M_p_ext(N+1) * p_undrained_poroelas_N; % ← was 0
    
    p_ext = thomas_algorithm(p_a, p_b, p_c, p_d, N);
    p     = p_ext(2:N-1);     % extract n physical-node values
    p = p'; T = T';

end
