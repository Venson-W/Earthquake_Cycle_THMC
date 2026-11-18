function [T, p] = implicit_p_T(Tn, pn, gamma, tauf, dt, T_boundarynorth, T_boundarysouth, p_boundarynorth, p_boundarysouth, M_alpha_hy, alpha_th, rho_c, Lambda, dy, n_pl_term, use_plastic_term, decomposition,rect_m_frac,Pr, Er, boundary)
    % Dirichlet keeps the production interior equations unchanged.
    if nargin >= 21 && strcmpi(boundary.type,'robin')
        [T,p]=legacy_robin_p_T(Tn,pn,gamma,tauf,dt, ...
            T_boundarynorth,T_boundarysouth,p_boundarynorth,p_boundarysouth, ...
            M_alpha_hy,alpha_th,rho_c,Lambda,dy,n_pl_term,use_plastic_term, ...
            decomposition,rect_m_frac,Pr,Er,boundary);
        return
    end
    % Get dimension - this is truly 1D, so Tn is a column vector
    ny = length(Tn);
    
    % Make sure we have at least 3 nodes (two boundaries + at least one interior)
    if ny < 3
        error('Grid must have at least 3 nodes for FEM solution');
    end
    
    % Initialize output arrays (must start with previous values)
    T = Tn; % Start with previous time step values
    p = pn; % Start with previous time step values
    
    % Pre-allocate arrays for tridiagonal system
    n_interior = ny-2;
    
    diag_a = zeros(n_interior, 1);
    diag_b = zeros(n_interior, 1);
    diag_c = zeros(n_interior, 1);
    diag_d = zeros(n_interior, 1);
    
    % Create system for all interior nodes (2 to ny-1)
    for i = 2:(ny-1)
        % Directly implement the FDM formula using implicit scheme
        % Original FDM: Te = tauf*gamma/rho_c + alpha_th*(T_{i+1}-2*T_i+T_{i-1})/dy^2
        % Implicit form: T_i - dt*alpha_th*(T_{i+1}-2*T_i+T_{i-1})/dy^2 = Tn_i + dt*tauf*gamma/rho_c
        
        % Set up a tridiagonal system
        a_diag = -dt * alpha_th / dy^2;  % Coefficient for T_{i-1}
        b_diag = 1 + 2 * dt * alpha_th / dy^2;  % Coefficient for T_i
        c_diag = -dt * alpha_th / dy^2;  % Coefficient for T_{i+1}
        
        % Right-hand side: Original value + source term
        % gamma is size ny-2, so we need to map i (2 to ny-1) to gamma index (1 to ny-2)
        gamma_idx = i - 1;
        d = Tn(i) + dt * (tauf * gamma(gamma_idx) / rho_c - rect_m_frac.*Er.*decomposition(gamma_idx));
        
        % For boundary nodes, include boundary values
        if i == 2
            d = d - a_diag * T_boundarysouth;  % Include south boundary
        end
        if i == ny-1
            d = d - c_diag * T_boundarynorth;  % Include north boundary
        end
        
        % Store coefficients in arrays for tridiagonal solver
        diag_a(i-1) = a_diag;
        diag_b(i-1) = b_diag;
        diag_c(i-1) = c_diag;
        diag_d(i-1) = d;
    end
    
    % Solve tridiagonal system for interior points
    % Using Thomas algorithm
    x = thomas_algorithm(diag_a, diag_b, diag_c, diag_d, n_interior);
    
    % Update temperature with solution
    for i = 2:(ny-1)
        T(i) = x(i-1);
    end
    
    % Apply boundary conditions
    T(1) = T_boundarysouth;
    T(ny) = T_boundarynorth;
    
    % ===== PRESSURE EQUATION =====
    % Calculate delta_T as in original FDM
    delta_T = (T - Tn) / dt;
    
    % Pre-allocate arrays for pressure tridiagonal system
    p_diag_a = zeros(n_interior, 1);
    p_diag_b = zeros(n_interior, 1);
    p_diag_c = zeros(n_interior, 1);
    p_diag_d = zeros(n_interior, 1);
    
    % Create system for all interior nodes (2 to ny-1)
    for i = 2:(ny-1)
        % Original FDM: pe = Lambda*delta_T + (M_alpha_hy_{i+1}*(p_{i+1}-p_i) - M_alpha_hy_i*(p_i-p_{i-1}))/dy^2
        % Implicit form: p_i - dt*(M_alpha_hy_{i+1}*(p_{i+1}-p_i) - M_alpha_hy_i*(p_i-p_{i-1}))/dy^2 = pn_i + dt*Lambda*delta_T
        
        % Coefficients for implicit form
        % M_alpha_hy has size ny-1
        m_idx = i - 1;  % Map to M_alpha_hy index (i-1)
        m_idx_plus = i;  % Map to M_alpha_hy index (i)
        
        a_coef = -dt * M_alpha_hy(m_idx) / dy^2;  % Coefficient for p_{i-1}
        c_coef = -dt * M_alpha_hy(m_idx_plus) / dy^2;  % Coefficient for p_{i+1}
        b_coef = 1 - a_coef - c_coef;  % Coefficient for p_i
        
        % Right-hand side: Original value + coupling term + optional plastic deformation term
        % Need to account for gamma and gamman size (ny-2)
        gamma_idx = i - 1;
        
        % Initialize the plastic deformation term
        plastic_term = 0;
        
        % Calculate the plastic deformation term from the updated FDM equation if switch is on
        if use_plastic_term
            plastic_term = n_pl_term(gamma_idx);
        end
        
        d_val = pn(i) + dt * (Lambda * delta_T(i) + plastic_term + rect_m_frac*Pr*decomposition(i));
        
        % For boundary nodes, include boundary values
        if i == 2
            d_val = d_val - a_coef * p_boundarysouth;  % Include south boundary
        end
        if i == ny-1
            d_val = d_val - c_coef * p_boundarynorth;  % Include north boundary
        end
        
        % Store coefficients in arrays for tridiagonal solver
        p_diag_a(i-1) = a_coef;
        p_diag_b(i-1) = b_coef;
        p_diag_c(i-1) = c_coef;
        p_diag_d(i-1) = d_val;
    end
    
    % Solve tridiagonal system for pressure
    p_x = thomas_algorithm(p_diag_a, p_diag_b, p_diag_c, p_diag_d, n_interior);
    
    % Update pressure with solution
    for i = 2:(ny-1)
        p(i) = p_x(i-1);
    end
    
    % Apply boundary conditions
    p(1) = p_boundarysouth;
    p(ny) = p_boundarynorth;
