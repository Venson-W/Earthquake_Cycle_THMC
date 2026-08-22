%% ========================================================================
%  Slider_ms.m — Coupled THMC Earthquake Cycle Simulator
% =========================================================================
%  A 1-D spring-slider model with a finite-width gouge layer embedded in
%  porous bulk material.  The code couples:
%    - Thermal  (frictional heating, thermal diffusion, thermal decomposition)
%    - Hydraulic (pore-pressure diffusion, permeability evolution)
%    - Mechanical (rate-and-state friction with a cross-scale RSF formulation)
%    - Chemical  (reaction kinetics for mineral decomposition)
%
%  Numerical method : Staggered Finite Difference Method (FDM) for the
%                     gouge layer, with Robin (diffusive) boundary conditions
%                     coupling to the semi-infinite porous bulk.
%
%  Developed by Yu-Han Wang 
%
%  Copyright 2025 Yu-Han Wang
%
%  Permission is hereby granted, free of charge, to any person obtaining a
%  copy of this software and associated documentation files (the "Software"),
%  to deal in the Software without restriction, including without limitation
%  the rights to use, copy, modify, merge, publish, distribute, sublicense,
%  and/or sell copies of the Software, and to permit persons to whom the
%  Software is furnished to do so, subject to the following conditions:
%
%  The above copyright notice and this permission notice shall be included
%  in all copies or substantial portions of the Software.
%
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
%  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
%  OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
%  ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%  OTHER DEALINGS IN THE SOFTWARE.
% =========================================================================

%% Initialization
clear; clc; close all; warning off;
addpath(fullfile(pwd, 'functions_ms'));
addpath(fullfile(pwd, 'postprocessing'));
addpath(fullfile(pwd, 'mat'));

% ---- Checkpoint / restart configuration -----------------------------------
%  Set RESTART = true to resume from a previously saved checkpoint file.
%  The checkpoint is written every CHECKPOINT_INTERVAL export steps and
%  upon normal termination or blow-up, so very little progress is lost.
RESTART              = false;
CHECKPOINT_FILE      = fullfile('..', 'post', 'checkpoint.mat');
CHECKPOINT_INTERVAL  = 500*25*100;   % save checkpoint every N export steps

% ---- Output directory -----------------------------------------------------
postFolder = fullfile('..', 'post');
if ~RESTART
    if exist(postFolder, 'dir')
        delete(fullfile(postFolder, '*'));
        disp('Post folder cleared.');
    else
        mkdir(postFolder);
        disp('Post folder created.');
    end
else
    if ~exist(postFolder, 'dir')
        error('Post folder does not exist — cannot restart.');
    end
end

% Reproducibility
rng(334500);

%% ========================================================================
%  Material Properties
% =========================================================================

% --- Elasticity -----------------------------------------------------------
G       = 30e9;           % Shear modulus [Pa]
density = 2700;           % Rock density  [kg/m^3]
cs      = sqrt(G / density);   % Shear-wave speed [m/s]
eta     = G / (2 * cs);        % Radiation-damping coefficient [Pa·s/m]
t_pull  = 0;                   % Time at which remote loading begins [s]
Vp      = @(t) (t > t_pull) * 1e-6;   % Plate loading velocity [m/s]
mu      = 1e-3;                % Fluid dynamic viscosity [Pa·s]

% --- Depth-dependent THMC parameters --------------------------------------
depth           = 7;      % Centroidal depth [km]
contrast_ratio  = 1;      % Permeability contrast between near-fault and bulk

% --- Physics switches (1 = on, 0 = off) -----------------------------------
permeability_evol    = 0;  % Porosity-dependent permeability evolution
decomposition_switch = 0;  % Thermal decomposition (reaction kinetics)
dilatancy            = 1;  % Dilatancy-induced pore-pressure changes

% --- Reaction kinetics parameters -----------------------------------------
rect_m_frac = 0.1;
A           = 10^15.47;   % Pre-exponential factor [1/s]
Q           = 319e3;      % Activation energy [J/mol]
R           = 8.3145;     % Universal gas constant [J/(K·mol)]
Pr          = 7.42e9;
rhof        = 418;
Tr          = 960;
Er          = 3060;
v           = 0.24;       % Poisson's ratio

% --- Depth-dependent parameters -------------------------------------------
if depth == 1
    Lambda   = 0.12e6;     % Thermal pressurization coefficient [Pa/K]
    rho_c    = 2.7e6;      % Volumetric heat capacity [J/(m^3·K)]
    alpha_th = 0.70e-6;    % Thermal diffusivity [m^2/s]
    alpha_hy = 7.15e-6;    % Hydraulic diffusivity [m^2/s]
    eps_pl   = 1.7e-4;     % Dilatancy coefficient [-]
    beta     = 4.39e-10;   % Pore-fluid compressibility [1/Pa]
    beta_T   = 1.75e-4;
    sigma_n  = 70e6;       % Normal stress [Pa]
    Ta       = 250;        % Ambient temperature [K]
    pa       = 10e6;       % Ambient pore pressure [Pa]
    k_permeability  = alpha_hy * beta * mu;
    k_conductivity_T = alpha_th * beta_T * mu;

elseif depth == 7
    Lambda   = 0.10e6;
    rho_c    = 2.7e6;
    alpha_th = 0.54e-6;
    alpha_hy = 6.71e-6;
    eps_pl   = 1.7e-4;
    beta     = 2.97e-10;
    beta_T   = 1.75e-4;
    sigma_n  = 140e6;
    pa       = 70e6;
    Ta       = 170;
    k_permeability  = alpha_hy * beta * mu;
    k_conductivity_T = alpha_th * beta_T * mu;
end

sigma_n0 = sigma_n;

% --- Initial porosity and derived hydro-mechanical parameters -------------
phi_0  = 0.05;
phi_e  = phi_0;
phi_pd = 0;
phi_pc = 0;
phi_p  = 0;
Z      = phi_0 / (1 - phi_0);

[beta_n, Lambda_n, beta_f, Lambda_f, eta_f, k_p] = ...
    invert_hydro_params(phi_0, Lambda, beta, alpha_hy);

beta_nu = beta_n / 2;
beta_fu = beta_f / 2;
beta_g  = 0.02e-9;
beta_gu = 0.011e-9;

% Maximum temperature under undrained-adiabatic thermal pressurization
Tmax = (1 + sqrt(alpha_hy / alpha_th)) * ((sigma_n - pa) / Lambda);
fprintf('Tmax (undrained-adiabatic) = %.1f K\n', Tmax);

%% ========================================================================
%  Rate-and-State Friction (RSF) Parameters
% =========================================================================
a    = 0.01;              % Direct-effect parameter
b    = a + 5e-3;          % Evolution-effect parameter  (velocity-weakening: b > a)
Drs  = 1e-3;              % State-evolution slip distance [m]
f0   = 0.6;               % Reference friction coefficient

%% ========================================================================
%  Spatial Discretization
% =========================================================================

% Characteristic RSF diffusion width (thermal-pressurization length scale)
W_rsf = pi^2 * a / f0^2 * rho_c / Lambda * (alpha_hy + alpha_th);

h   = 0.0005;             % Gouge-layer half-width [m]
ny  = 300;                % Number of grid points across the gouge layer
nx  = 1;                  % Spring-slider: single spatial element

yy          = linspace(-h/2, h/2, ny);
yy_gammadot = yy;
dy          = h / (ny - 1);

%% ========================================================================
%  Initial Conditions
% =========================================================================

% Background shear stress (at steady-state friction)
taub0 = f0 * (sigma_n - pa);
tauf  = taub0;

% Initial field variables
p  = zeros(1, ny) + pa * randn(1, ny) * 0.001;   % Small random perturbation
T  = zeros(1, ny);
xi = zeros(1, ny);          % Reaction progress (0 = unreacted, 1 = fully reacted)
dx = 1e9;

% Characteristic slip for pressurization under undrained-adiabatic conditions
deltac = rho_c * h / (f0 * Lambda);

% --- Spring stiffness and stability check ---------------------------------
k_crit = abs((sigma_n - pa) * (b - a) / Drs);   % Critical stiffness
R_k    = 0.5;                                     % Stiffness ratio k/k_crit
k      = k_crit * R_k;
L      = G / k;
K      = 2 * G * (1 - v) / (1 - 2*v);            % Bulk modulus
kn     = K / L;                                    % Normal stiffness

if R_k >= 1
    disp('Unconditionally stable (R_k >= 1).');
else
    fprintf('Conditionally unstable: R_k = %.3f\n', R_k);
end

% Block initial position and velocity
x_0         = 0;
x           = x_0;
V0          = 1e-6;           % Reference slip velocity [m/s]
V_slider0   = V0;
gammadot0   = V0 / h;         % Reference strain rate [1/s]
gamma       = 0;              % Accumulated shear strain
V           = V_slider0;
gammadot    = (V / h) * ones(1, ny);

% State variable (ageing law)
theta0  = Drs / V0;
theta_r = theta0;
theta   = theta0 * ones(1, ny);

n_pl_term = zeros(1, ny);

% --- Hydraulic diffusivity mobility matrix --------------------------------
%  Harmonic averaging is used on staggered grid; near-boundary nodes can
%  have a permeability contrast relative to the interior.
M_alpha_hy = alpha_hy * ones(1, ny+1);
M_alpha_hy(:, 1:2)       = alpha_hy / contrast_ratio;
M_alpha_hy(:, end-3:end) = alpha_hy / contrast_ratio;

M_kappa_y = ones(nx+1, ny+1) * alpha_hy * beta;
M_kappa_y(:, 1:2)       = M_kappa_y(:, 1:2)       / contrast_ratio;
M_kappa_y(:, end-3:end) = M_kappa_y(:, end-3:end) / contrast_ratio;

M_alpha_hy0 = M_alpha_hy;     % Store initial mobility for toggle

% Boundary conditions (Robin / diffusive — values at infinity)
[T_boundarynorth, T_boundarysouth, ...
 p_boundarynorth, p_boundarysouth] = deal(0, 0, 0, 0);

%% ========================================================================
%  Bulk Diffusion Domain (geometrically stretched grid)
% =========================================================================
refine = 1;
n_bulk = 35;
n_bulk = refine * n_bulk + 1;
W_bulk = 100;   % Total width of the bulk domain [m]

% Solve for the geometric stretch ratio r so that the first grid spacing
% equals dy (matching the gouge grid) and the total extent equals W_bulk.
f_stretch = @(r) sum(dy * r.^(0:n_bulk-2)) - W_bulk;
r_bulk    = fzero(f_stretch, 1.01);
dy_bulk   = dy * r_bulk.^(0:n_bulk-2);
y_bulk    = 1 + [0, cumsum(dy_bulk)]';

%% ========================================================================
%  Time-Stepping Configuration
% =========================================================================
dt              = 1e-7;                 % Initial time step [s]
secs_per_year   = 365 * 24 * 3600;
secs_per_hour   = 3600;
current_time    = 0;
current_time_year = 0;
count           = 0;

ending_time     = 3600 * 24 * 15;      % Total simulation time [s]
dt_limit        = ending_time / 1e5;    % Maximum allowed time step [s]

% --- Iterative solver parameters ------------------------------------------
max_iter   = 50;          % Maximum sub-iterations per time step
miniter    = 1;           % Minimum sub-iterations before accepting
frac       = 0.5 / (2^4);
maxtolvio  = 0;
maxdtfac   = 1.05;       % Maximum time-step growth factor

% --- Colour palette (for optional live plotting) --------------------------
alphaVal = 0.25;
c1 = [0.0000, 0.4470, 0.7410, 1.00];      % Blue
c2 = [0.8500, 0.3250, 0.0980, 1.00];      % Orange
c3 = [0.4660, 0.6740, 0.1880, alphaVal];  % Green
c4 = [0.4940, 0.1840, 0.5560, alphaVal];  % Purple
m_to_mm = 1e3;

% --- Store previous-step values for the predictor-corrector ---------------
gammadot_p = gammadot;
V_p        = V;
theta_p    = theta;
dtp        = dt;
x_p        = x;
p_p        = p;

%% ========================================================================
%  Pre-allocate History Arrays
% =========================================================================
max_steps = 3e6;

timev               = zeros(max_steps, 1);
history_V           = zeros(max_steps, 1);
history_x           = zeros(max_steps, 1);
history_tauf        = zeros(max_steps, 1);
history_state       = zeros(max_steps, length(theta));
history_xi          = zeros(max_steps, length(xi));
history_gammadot    = zeros(max_steps, length(gammadot));
history_gamma       = zeros(max_steps, length(gammadot));
history_T           = zeros(max_steps, length(T));
history_p           = zeros(max_steps, length(p));
history_maxgammadot = zeros(max_steps, 1);
history_p_max       = zeros(max_steps, 1);
history_p_avg       = zeros(max_steps, 1);
history_p_r         = zeros(max_steps, 1);
history_f_r         = zeros(max_steps, 1);
timev_secs          = zeros(max_steps, 1);
history_i_export    = zeros(max_steps, 1);
history_phi         = zeros(max_steps, length(gammadot));
history_phi_e       = zeros(max_steps, length(gammadot));
history_phi_p       = zeros(max_steps, length(gammadot));
history_alphahy     = zeros(max_steps, length(gammadot) + 1);
history_theta_r     = zeros(max_steps, 1);

i_export = 0;

%% ========================================================================
%  Normalisation Scales
% =========================================================================
gammadotc = 2 * cs * taub0 / (G * h);               % Dynamic strain-rate scale
fc        = f0 + (a - b) * log(gammadotc / gammadot0);  % Dynamic friction scale
deltac    = rho_c * h / Lambda / fc;                 % Slip scale for pressurization
tc        = 2 * cs * taub0 / (G * deltac);           % Dynamic time scale
gammac    = Drs / h;                                 % Characteristic strain

phi = phi_0 * ones(1, ny);

%% ========================================================================
%  Restore Checkpoint (if restarting)
% =========================================================================
if RESTART
    if ~isfile(CHECKPOINT_FILE)
        error('Checkpoint file not found: %s', CHECKPOINT_FILE);
    end
    fprintf('Loading checkpoint from %s ...\n', CHECKPOINT_FILE);
    ckpt = load(CHECKPOINT_FILE);

    % --- Recover simulation state -----------------------------------------
    current_time      = ckpt.current_time;
    current_time_year = ckpt.current_time_year;
    count             = ckpt.count;
    dt                = ckpt.dt;
    dtp               = ckpt.dtp;

    % Field variables
    V           = ckpt.V;
    V_p         = ckpt.V_p;
    x           = ckpt.x;
    x_p         = ckpt.x_p;
    p           = ckpt.p;
    p_p         = ckpt.p_p;
    T           = ckpt.T;
    theta       = ckpt.theta;
    theta_p     = ckpt.theta_p;
    theta_r     = ckpt.theta_r;
    xi          = ckpt.xi;
    tauf        = ckpt.tauf;
    gamma       = ckpt.gamma;
    gammadot    = ckpt.gammadot;
    gammadot_p  = ckpt.gammadot_p;
    phi         = ckpt.phi;
    phi_e       = ckpt.phi_e;
    phi_p       = ckpt.phi_p;
    phi_pc      = ckpt.phi_pc;
    phi_pd      = ckpt.phi_pd;
    sigma_n     = ckpt.sigma_n;
    M_alpha_hy  = ckpt.M_alpha_hy;
    k_p         = ckpt.k_p;
    n_pl_term   = ckpt.n_pl_term;
    maxtolvio   = ckpt.maxtolvio;

    % --- Recover history arrays -------------------------------------------
    i_export            = ckpt.i_export;
    timev_secs(1:ckpt.count+1) = ckpt.timev_secs(1:ckpt.count+1);

    timev(1:i_export)               = ckpt.timev(1:i_export);
    history_V(1:i_export)           = ckpt.history_V(1:i_export);
    history_x(1:i_export)           = ckpt.history_x(1:i_export);
    history_tauf(1:i_export)        = ckpt.history_tauf(1:i_export);
    history_state(1:i_export,:)     = ckpt.history_state(1:i_export,:);
    history_p(1:i_export,:)         = ckpt.history_p(1:i_export,:);
    history_gammadot(1:i_export,:)  = ckpt.history_gammadot(1:i_export,:);
    history_gamma(1:i_export+1,:)   = ckpt.history_gamma(1:i_export+1,:);
    history_maxgammadot(1:i_export) = ckpt.history_maxgammadot(1:i_export);
    history_T(1:i_export,:)         = ckpt.history_T(1:i_export,:);
    history_xi(1:i_export,:)        = ckpt.history_xi(1:i_export,:);
    history_p_max(1:i_export)       = ckpt.history_p_max(1:i_export);
    history_p_avg(1:i_export)       = ckpt.history_p_avg(1:i_export);
    history_phi(1:i_export,:)       = ckpt.history_phi(1:i_export,:);
    history_phi_e(1:i_export,:)     = ckpt.history_phi_e(1:i_export,:);
    history_phi_p(1:i_export,:)     = ckpt.history_phi_p(1:i_export,:);
    history_alphahy(1:i_export,:)   = ckpt.history_alphahy(1:i_export,:);
    history_p_r(1:i_export)         = ckpt.history_p_r(1:i_export);
    history_f_r(1:i_export)         = ckpt.history_f_r(1:i_export);
    history_i_export(1:i_export)    = ckpt.history_i_export(1:i_export);
    history_theta_r(1:i_export)     = ckpt.history_theta_r(1:i_export);

    clear ckpt;
    fprintf('  Resumed at t = %.6e s  (count = %d, i_export = %d)\n', ...
            current_time, count, i_export);
end

%% ========================================================================
%  Main Time-Stepping Loop
% =========================================================================
tic;

while current_time < ending_time

    % ------------------------------------------------------------------
    %  Predictor–corrector iteration control
    % ------------------------------------------------------------------
    tryagaincount = 0;
    enter         = 1;

    % Explicit predictor (linear extrapolation from two previous steps)
    V_g        = V        + (V        - V_p)        / dtp * dt;
    V_mg       = 0.5 * (V + V_g);
    gammadot_g = gammadot + (gammadot - gammadot_p) / dtp * dt;
    x_g        = x        + (x        - x_p)        / dtp * dt;
    p_g        = p        + (p        - p_p)        / dtp * dt;

    % Store current-step values (will become "previous" after acceptance)
    V_p              = V;
    gammadot_p       = gammadot;
    theta_p          = theta;
    x_p              = x;
    tauf_p           = tauf;
    p_p              = p;
    T_p              = T;
    phi_ep           = phi_e;
    phi_pcp          = phi_pc;
    phi_pdp          = phi_pd;
    phi_pp           = phi_p;
    xi_p             = xi;
    current_time_p   = current_time;
    count            = count + 1;

    % ------------------------------------------------------------------
    %  Sub-iteration loop (adaptive time-step halving)
    % ------------------------------------------------------------------
    while enter || tryagaincount

        switch tryagaincount
            case 0
                % Use the explicit predictor as the initial guess.
            case 1
                % First corrector: refine without halving dt.
                p_g    = p;
                T_g    = T;
                x_g    = x;
                V_g    = V;
                V_mg   = 0.5 * (V_g + V_p);
            otherwise
                % Subsequent corrections: average old and new.
                p_g    = 0.5 * (p + p_p);
                T_g    = 0.5 * (T + T_p);
                x_g    = 0.5 * (x + x_p);
                V_g    = 0.5 * (V_p + V);
                V_mg   = 0.5 * (V_p + V_g);
        end

        enter = 0;

        % --- Mid-step estimates (Crank–Nicolson style) ----------------
        V_m        = 0.5 * (V_g + V_p);
        gammadot_m = 0.5 * (gammadot_p + gammadot_g);

        % --- State-variable update (ageing law) -----------------------
        theta_predictor = state_theta(theta_p, V_m, dt, Drs);

        % --- Slip predictor -------------------------------------------
        x_predictor = x_p + V_m * dt;

        % --- Solve for fault shear strength ---------------------------
        [tauf, nb_it] = solvestrength( ...
            tauf, sigma_n, p + pa, theta, ...
            k * (Vp(current_time + dt) * (current_time + dt) - x_predictor) + taub0, ...
            gammadot0, V0, a, f0, b, Drs, eta, yy_gammadot);

        % --- Thermal decomposition (Arrhenius kinetics) ---------------
        decomposition = (1 - xi) .* A .* exp(-Q ./ (R .* (T_p + Ta))) ...
                        .* decomposition_switch;

        % --- Dilatancy source term for pore-pressure equation ---------
        n_pl_term = gammadot_p ./ gammac ...
                    .* (phi_pdp - eps_pl .* log(gammadot_p ./ gammadot0)) ...
                    * (1 / beta) .* dilatancy;

        % --- Implicit solve: coupled pore-pressure & temperature ------
        [T, p] = implicit_p_T( ...
            T_p, p_p, ...
            gammadot, tauf, dt, ...
            T_boundarynorth, T_boundarysouth, ...
            p_boundarynorth, p_boundarysouth, ...
            M_alpha_hy, alpha_th, rho_c, Lambda, dy, ...
            n_pl_term, dilatancy, decomposition, ...
            rect_m_frac, Pr, Er, alpha_hy, alpha_th);

        % --- Porosity evolution (elastic + plastic) -------------------
        phi_e = phi_ep + phi_0 .* (beta_n * (p - p_p) + Lambda_n * (T - T_p));
        phi   = phi_e + phi_p;

        % --- Permeability update (cubic porosity law) -----------------
        k_pnew = k_p .* (phi / phi_0).^3;

        % --- Update mobility matrix with evolved permeability ---------
        M_alpha_hy = [k_pnew(1), ...
                      2 ./ (1./k_pnew(1:end-1) + 1./k_pnew(2:end)), ...
                      k_pnew(end)] ...
                     ./ beta ./ eta_f .* permeability_evol ...
                     + M_alpha_hy0 .* (1 - permeability_evol);

        % --- Update reaction progress ---------------------------------
        if decomposition_switch == 1
            [xi, kDt] = reaction_progress(xi_p, T_p + Ta, T, dt, A, Q, R);
        end

        % --- Strain rate from cross-scale RSF constitutive law --------
        gammadot = shear_strain( ...
            tauf, theta_predictor, gammadot0, V0, f0, b, a, Drs, ...
            sigma_n, p + pa);

        % --- Integrate slip velocity from strain-rate profile ---------
        V = slip_velocity(gammadot, yy_gammadot);

        % --- Corrector: mid-point averaging and slip update -----------
        V_m = 0.5 * (V + V_p);
        x   = x_p + V_m * dt;

        % --- Poroelastic normal-stress coupling -----------------------
        Dy      = trapz(yy, ...
                    (Z .* beta_n - beta_g) .* ...
                    (p - (Z .* beta_nu + beta_gu) ./ (Z .* beta_n - beta_g) .* sigma_n) ...
                    + phi_p ./ (1 - phi_0), 2);
        sigma_n = sigma_n0 + kn * Dy;

        % --- State update (ageing law, corrected mid-point) -----------
        theta = state_theta(theta_p, V_m, dt, Drs);

        % --- Plastic porosity: chemical compaction --------------------
        phi_pc = phi_pcp ...
                 - (1 - xi) .* A .* exp(-Q ./ (R .* (T + Ta))) ...
                 .* beta .* decomposition_switch;

        % --- Plastic porosity: shear-enhanced dilatancy ---------------
        phi_pd = (phi_pdp / dt ...
                  + gammadot ./ gammac .* (eps_pl .* log(gammadot ./ gammadot0))) ...
                 ./ (1/dt + gammadot ./ gammac) .* dilatancy;

        phi_p = phi_pd + phi_pc;

        shear_stress = tauf;

        % --- Convergence check ----------------------------------------
        [accuracyabs, ~] = max(abs(p - p_g) ./ (a .* sigma_n));

        if (tryagaincount < miniter ...
                || max(abs(p - p_g) ./ (a .* sigma_n)) > 0.1 * frac ...
                || norm(p - p_g) / (norm(p) + 1) > 0.1 * frac) ...
                && tryagaincount < max_iter

            tryagain      = 1;
            tryagaincount = tryagaincount + 1;

            if tryagaincount == 1
                % First iteration: refine prediction without halving dt
            else
                dt = dt / 2;
            end
        else
            tryagain            = 0;
            tryagaincount_hist  = tryagaincount;
            tryagaincount       = 0;

            if accuracyabs > maxtolvio
                maxtolvio     = accuracyabs;
                maxtolviotime = current_time;
            end
        end

        % --- Accept time step and compute next dt ---------------------
        if tryagain == 0
            current_time      = current_time_p + dt;
            current_time_year = current_time / secs_per_year;

            C1        = 1/16;
            dt_theta1 = C1 * Drs / V;
            dtp       = dt;
            dt        = min([dt_theta1, dt_limit]);

            if dt > maxdtfac * dtp
                dt = maxdtfac * dtp;
            end
        end

    end  % sub-iteration loop

    timev_secs(count + 1) = current_time;

    %% ------------------------------------------------------------------
    %  Compute Representative (Fault-Averaged) Quantities
    % -------------------------------------------------------------------
    p_r = sigma_n - shear_stress / (f0 + a * log(V_m / V0) ...
          + b * log(theta_r / theta0)) - pa;
    f_r = shear_stress / (sigma_n - p_r - pa);

    gamma = gamma + gammadot * dt;

    %% ------------------------------------------------------------------
    %  Save History (every 25 steps or after difficult convergence)
    % -------------------------------------------------------------------
    if (mod(count, 25) == 0) || (tryagaincount_hist > 2)

        i_export = i_export + 1;

        timev(i_export)               = real(current_time_year);
        history_V(i_export)           = real(V);
        history_x(i_export)           = real(x);
        history_tauf(i_export)        = real(shear_stress);
        history_state(i_export,:)     = real(theta);
        history_p(i_export,:)         = real(p);
        history_gammadot(i_export,:)  = real(gammadot);
        gammadot                      = real(gammadot);
        history_gamma(i_export+1,:)   = gamma(:);
        history_maxgammadot(i_export) = max(gammadot);
        history_T(i_export,:)         = T;
        history_xi(i_export,:)        = xi;
        history_p_max(i_export)       = max(p);
        history_p_avg(i_export)       = mean(p);
        history_phi(i_export,:)       = phi;
        history_phi_e(i_export,:)     = phi_e;
        history_phi_p(i_export,:)     = phi_p;
        history_alphahy(i_export,:)   = M_alpha_hy;
        history_p_r(i_export)         = p_r;
        history_f_r(i_export)         = f_r;
        history_i_export(i_export)    = i_export;

        % --- Periodic checkpoint save ---------------------------------
        if mod(i_export, CHECKPOINT_INTERVAL) == 0
            save_checkpoint(CHECKPOINT_FILE, ...
                current_time, current_time_year, count, dt, dtp, ...
                V, V_p, x, x_p, p, p_p, T, theta, theta_p, theta_r, ...
                xi, tauf, gamma, gammadot, gammadot_p, ...
                phi, phi_e, phi_p, phi_pc, phi_pd, sigma_n, ...
                M_alpha_hy, k_p, n_pl_term, maxtolvio, ...
                i_export, timev_secs, timev, ...
                history_V, history_x, history_tauf, history_state, ...
                history_p, history_gammadot, history_gamma, ...
                history_maxgammadot, history_T, history_xi, ...
                history_p_max, history_p_avg, history_phi, ...
                history_phi_e, history_phi_p, history_alphahy, ...
                history_p_r, history_f_r, history_i_export, ...
                history_theta_r);
            fprintf('  [Checkpoint] saved at t = %.6e s  (export #%d)\n', ...
                    current_time, i_export);
        end
    end

    %% ------------------------------------------------------------------
    %  Blow-up Detection
    % -------------------------------------------------------------------
    if max(real(V)) > 100 || any(isnan(V)) || any(isinf(V))
        disp('Simulation terminated: blow-up or invalid values detected.');
        % Save a final checkpoint so the user can inspect state
        save_checkpoint(CHECKPOINT_FILE, ...
            current_time, current_time_year, count, dt, dtp, ...
            V, V_p, x, x_p, p, p_p, T, theta, theta_p, theta_r, ...
            xi, tauf, gamma, gammadot, gammadot_p, ...
            phi, phi_e, phi_p, phi_pc, phi_pd, sigma_n, ...
            M_alpha_hy, k_p, n_pl_term, maxtolvio, ...
            i_export, timev_secs, timev, ...
            history_V, history_x, history_tauf, history_state, ...
            history_p, history_gammadot, history_gamma, ...
            history_maxgammadot, history_T, history_xi, ...
            history_p_max, history_p_avg, history_phi, ...
            history_phi_e, history_phi_p, history_alphahy, ...
            history_p_r, history_f_r, history_i_export, ...
            history_theta_r);
        break;
    end

end  % main time loop

%% ========================================================================
%  Post-Processing: Trim & Save
% =========================================================================
t_wallclock = toc;

% Trim history arrays to actual length
timev               = timev(1:i_export);
history_V           = history_V(1:i_export);
history_x           = history_x(1:i_export);
history_tauf        = history_tauf(1:i_export);
history_maxgammadot = history_maxgammadot(1:i_export);
run_time            = t_wallclock;
history_p_r         = history_p_r(1:i_export);
history_f_r         = history_f_r(1:i_export);
history_p_max       = history_p_max(1:i_export);
history_theta_r     = history_theta_r(1:i_export);

history_gamma       = history_gamma(1:i_export,:);
history_gammadot    = history_gammadot(1:i_export,:);
history_p           = history_p(1:i_export,:);
history_T           = history_T(1:i_export,:);
history_state       = history_state(1:i_export,:);
history_xi          = history_xi(1:i_export,:);
history_phi         = history_phi(1:i_export,:);
history_phi_e       = history_phi_e(1:i_export,:);
history_phi_p       = history_phi_p(1:i_export,:);
history_alphahy     = history_alphahy(1:i_export,:);

% --- Save results ---------------------------------------------------------
save(fullfile('..', 'post', 'fault_slip_data.mat'), ...
    'timev', 'history_V', 'history_f_r', 'history_p_max', ...
    'history_x', 'history_tauf', 'history_maxgammadot', 'run_time', ...
    'a', 'f0', 'Lambda', 'alpha_th', 'alpha_hy', 'rho_c', ...
    'ny', 'h', 'yy', 'gammadotc', 'deltac', 'fc', 'tc', '-v7');

save(fullfile('..', 'post', 'fault_slip_field.mat'), ...
    'history_gammadot', 'history_p', 'history_T', 'history_state', ...
    'history_xi', 'history_phi', 'history_phi_e', 'history_phi_p', ...
    'history_alphahy', '-v7');

% --- Final checkpoint (for potential future continuation) -----------------
save_checkpoint(CHECKPOINT_FILE, ...
    current_time, current_time_year, count, dt, dtp, ...
    V, V_p, x, x_p, p, p_p, T, theta, theta_p, theta_r, ...
    xi, tauf, gamma, gammadot, gammadot_p, ...
    phi, phi_e, phi_p, phi_pc, phi_pd, sigma_n, ...
    M_alpha_hy, k_p, n_pl_term, maxtolvio, ...
    i_export, timev_secs, timev, ...
    history_V, history_x, history_tauf, history_state, ...
    history_p, history_gammadot, history_gamma, ...
    history_maxgammadot, history_T, history_xi, ...
    history_p_max, history_p_avg, history_phi, ...
    history_phi_e, history_phi_p, history_alphahy, ...
    history_p_r, history_f_r, history_i_export, ...
    history_theta_r);

fprintf('Simulation completed after %.1f seconds (wall-clock).\n', t_wallclock);


%% ========================================================================
%  Local Function: save_checkpoint
% =========================================================================
function save_checkpoint(filepath, ...
    current_time, current_time_year, count, dt, dtp, ...
    V, V_p, x, x_p, p, p_p, T, theta, theta_p, theta_r, ...
    xi, tauf, gamma, gammadot, gammadot_p, ...
    phi, phi_e, phi_p, phi_pc, phi_pd, sigma_n, ...
    M_alpha_hy, k_p, n_pl_term, maxtolvio, ...
    i_export, timev_secs, timev, ...
    history_V, history_x, history_tauf, history_state, ...
    history_p, history_gammadot, history_gamma, ...
    history_maxgammadot, history_T, history_xi, ...
    history_p_max, history_p_avg, history_phi, ...
    history_phi_e, history_phi_p, history_alphahy, ...
    history_p_r, history_f_r, history_i_export, ...
    history_theta_r)
% SAVE_CHECKPOINT  Write all simulation state to a MAT-file for restart.
%   Saves to a temporary file first, then atomically renames to avoid
%   corruption if the process is killed during the write.

    tmpfile = [filepath, '.tmp'];
    save(tmpfile, ...
        'current_time', 'current_time_year', 'count', 'dt', 'dtp', ...
        'V', 'V_p', 'x', 'x_p', 'p', 'p_p', 'T', ...
        'theta', 'theta_p', 'theta_r', ...
        'xi', 'tauf', 'gamma', 'gammadot', 'gammadot_p', ...
        'phi', 'phi_e', 'phi_p', 'phi_pc', 'phi_pd', 'sigma_n', ...
        'M_alpha_hy', 'k_p', 'n_pl_term', 'maxtolvio', ...
        'i_export', 'timev_secs', 'timev', ...
        'history_V', 'history_x', 'history_tauf', 'history_state', ...
        'history_p', 'history_gammadot', 'history_gamma', ...
        'history_maxgammadot', 'history_T', 'history_xi', ...
        'history_p_max', 'history_p_avg', 'history_phi', ...
        'history_phi_e', 'history_phi_p', 'history_alphahy', ...
        'history_p_r', 'history_f_r', 'history_i_export', ...
        'history_theta_r', ...
        '-v7');
    movefile(tmpfile, filepath);
end
