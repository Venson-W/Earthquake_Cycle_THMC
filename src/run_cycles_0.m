function result=run_cycles_0(options)
%% README
% Current manuscript case cycles_0. Software citation: Yuhan Wang.
% ORCID: https://orcid.org/0000-0003-0490-5977

% Please cite this software as specified in CITATION.cff.
% Original development by Yuhan Wang and Elias Rafn Heimisson.

%Copyright 2025 Yuhan Wang

%Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Problem setup
% 1d spring slider with finite width and diffusive boundaries, representing the gouge layer embedded in porous bulks. The gouge layer allows the fluid to diffuse in a any porous materials insofar as permeabilities are concerned. The bulk is considered homogeneous and isotropic in permeability. 
% Two numerical schemes are combined herein: 
% (i) Staggered Finite Difference Method, for the gouge layer; 
% (ii) Spectral Boundary Ingegral Method, for the diffusive bulk.

% clean history
% Preserve the caller workspace, figures and warning state.
inelast_options=options;
helper_path=fileparts(mfilename('fullpath'));
previous_path=path;cleanup_path=onCleanup(@()path(previous_path));
addpath(helper_path);

% Define the path to the post folder
postFolder = char(options.output_directory);

assert(~isfolder(postFolder) && ~isfile(postFolder),'Use a fresh output folder.');
mkdir(postFolder);

% Random seed
rng(334500); % Set random seed based on current time

%% Material properties
% Elasticity
G = 30e9;
density = 2700;             
cs = sqrt(G / density); 
eta = G / (2 * cs);
L = 60 * 1000; % k = G / L; Spring stiffness
t_pull = 0;
Vp = @(t) (t > t_pull) * 1e-6;
mu = 1e-3;

% Thermal pressurization
depth = 7; % centroidal depth (km)
contrast_ratio = 1;

% Permeability evolution switch
permeability_evol = 0;
% Thermal decomposition switch
decomposition_switch = 0;
% Dilatancy switch
dilatancy = 0;

rect_m_frac = 0.5;
A = 10^15.47;
Q = 319e3;
R = 8.3145; % the gas constant [R] = J/K/mol
Pr = 7.42e9;
rhof = 418;
Tr = 960; 
Er = 3060;

if depth == 1
    Lambda = 0.12e6; % Thermal pressurization
    rho_c = 2.7e6; % 
    alpha_th = 0.70e-6;
    alpha_hy = 7.15e-6; 
    eps_pl = 1.7e-4; %
    beta = 4.39e-10; % Compressibility (1/Pa)  
    beta_T = 1.75e-4; % 3e-6
    sigma_n = 70e6;
    Ta = 250;
    pa = 10e6;
    k_permeability = alpha_hy * beta * mu;
    k_conductivity_T = alpha_th * beta_T * mu;
  
elseif depth == 7
    Lambda = 0.10e6; % Thermal pressurization
    rho_c = 2.7e6; % 
    alpha_th = 0.54e-6;
    alpha_hy = 6.71e-6; 
    eps_pl = 1.7e-4; % 
    beta = 2.97e-10; % Compressibility (1/Pa)  
    beta_T = 1.75e-4; % 3e-6
    sigma_n = 140e6;
    pa = 70e6;
    Ta = 170;
    k_permeability = alpha_hy * beta * mu;
    k_conductivity_T = alpha_th * beta_T * mu;
end

phi_0 = 0.05; phi_e = phi_0; phi_pd =0; phi_pc = 0; phi_p = 0;
[beta_n, Lambda_n, beta_f, Lambda_f, eta_f, k_p] = invert_hydro_params(phi_0, Lambda, beta, alpha_hy);

% Max T
Tmax = (1 + sqrt(alpha_hy/alpha_th))*((sigma_n - pa)/Lambda)

%% Numerical discretization
% Time
dt = 1e-7; % time step for the first count

secs_per_year = 365 * 24 * 60 * 60;      

secs_per_hour = 3600;

current_time = 0;

count = 0;

%% RSF parameters and initial stress condition 
a = 0.01;
b = a + 5e-3; %
Dc = 5e-3;
f0 = 0.6;


% Space grid
W_rsf = pi.^2 .* a ./ f0.^2 .*rho_c ./ Lambda .* (alpha_hy + alpha_th) ./1;
h = 0.0005;%W_rsf/10;

ny = 300; % Counts of the y-axis node

nx=1; % spring slider

yy = linspace(-h/2, h/2, ny);

yy_gammadot = yy;

dy = h/(ny-1);     % Y-axis grid size

% initial field

taub0 = f0 * (sigma_n - pa); % Background shear stress
tauf = taub0;

p = zeros(1,ny) + pa * randn(1, ny) * 0.001;
T = zeros(1,ny); 
xi = zeros(1,ny);
dx = 1e9;
%xi = 0.25 * ((a - b) ./ a).^2; %xi = 0.25 * (K * Dc ./ (a .* P) - (b - a) ./ a).^2 - K * Dc ./ (a .* p);

% Slip scale in which pressurization occurs under undrained, adiabatic conditions.
deltac = rho_c * h / (f0 * Lambda);

% Check spring stiffness k value and k_crit
k_crit = abs((sigma_n-pa) * (b - a) / Dc);
R_k = 0.5;
k = k_crit*R_k;

if (R_k > 1) || (R_k == 1)
    disp("unconditioanally stable");
else
    disp(strcat("Rk = ",num2str(R_k)));
end

% Initial condition
x_0 = 0; x = x_0; % initial location of block

V0 = 1e-6; 

V_slider0 = V0; % initial slip velocity of the block

gammadot0 = V0 ./ h; 

gamma = 0;

V = V_slider0;

gammadot = V./h .* ones(1,ny);
% Isolated steady-porosity reference; keep phi_pd(0)=0 and initial total porosity.
% This is the prescribed uniform initial rate BEFORE deliberate pressure noise.
gammadot_dilatancy_ref = (V_slider0/h) .* ones(size(gammadot));
inelast_reference_profile = 'uniform';


theta0 = Dc ./ V0;

theta_r = theta0; % representative state 

theta = theta0 .* ones(1,ny);

n_pl_term = zeros(1,ny);

% W_rsf = 6.9 .* (a-b) .* rho_c .* (sqrt(alpha_hy) + sqrt(alpha_th)) .^2 ./f0 ./Lambda ./V .* (f0+2*(a-b));

% Mobility matrix
M_alpha_hy = alpha_hy * ones(1, ny+1); M_alpha_hy(:,1:2) = alpha_hy./contrast_ratio; M_alpha_hy(:,end-3:end) = alpha_hy./contrast_ratio;   
M_kappa_y = ones(nx+1,ny+1).*alpha_hy*beta; M_kappa_y(:,1:2) = M_kappa_y(:,1:2)./contrast_ratio; M_kappa_y(:,end-3:end) = M_kappa_y(:,end-3:end)./contrast_ratio;
M_alpha_hy0 = M_alpha_hy;
% Boundary Condition, Dirichlet
[T_boundarynorth, T_boundarysouth, p_boundarynorth, p_boundarysouth] = deal(0, 0, 0, 0);

% Created for post-processing
current_time_year = 0;
%% Shear Zone Numerical Setup
max_iter = 50;

%% Define bulk nodes
% spatial discretization and diffusion parameters
refine = 1; % for the convenience of imposing finer grid in bulks
n_bulk = 35; 
n_bulk = refine*n_bulk + 1; 
W_bulk = 100; % width of bulks (m)
f = @(r) sum(dy * r.^(0:n_bulk-2)) - W_bulk; % the first grid size at one end of the bulk is "dy", the same as that of the fault gouge. Then a 
r_bulk = fzero(f, 1.01);   % solve for r (start guess slightly > 1)
dy_bulk = dy .* r_bulk.^(0:n_bulk-2); % grid sizes for the bulk
y_bulk = 1+[0, cumsum(dy_bulk)]'; % local coordinates, perpendicular to the fault
%y_bulk = linspace(0, n_bulk*dy, n_bulk)';

%%
ending_time = 3600*24*60;%(Dc/h)/gammadot0*1000;%secs_per_year * 1;
dt_limit = ending_time/1e5;
% Release options preserve the original defaults unless explicitly supplied.
if isfield(options,'end_seconds'),ending_time=options.end_seconds;end
if isfield(options,'dt_limit'),dt_limit=options.dt_limit;end
release_boundary=thmc_boundary_options(options,alpha_hy,alpha_th,rho_c,h);




% Define transparent colors with alpha = 0.75
alphaVal = 0.25;
c1 = [0, 0.4470, 0.7410, 1];   % Blue
c2 = [0.8500, 0.3250, 0.0980,1 ]; % Orange
c3 = [0.4660, 0.6740, 0.1880, alphaVal]; % Green
c4 = [0.4940, 0.1840, 0.5560, alphaVal]; % Purple
m_to_mm = 1e3;


% Initialize variables for the iterative process
miniter = 1;
frac = 0.5 ./ (2 ^ 4);
maxtolvio = 0;
inelast_audit=struct('pressure_tolerance',0.1*frac,'accepted_step_count',0, ...
    'nonconverged_steps',0,'invalid_state_steps',0, ...
    'maximum_absolute_normalized_pressure_residual',0,'maximum_relative_pressure_residual',0, ...
    'maximum_pressure_retries',0,'minimum_sigma_e',Inf,'minimum_gammadot',Inf,'minimum_theta',Inf);

maxdtfac = 1.05;

% Initialize previous values for predictor-corrector method
gammadot_p = gammadot;
V_p = V;
theta_p = theta;
dtp = dt;
x_p = x;
p_p = p;

% Valving
%eta_f = 1e-3;
%ah_min = 1e-19/eta_f/beta;
%ah_max = 1e-15/eta_f/beta;
%L_valving = 1;
%sigma_eff_star = 30e6;

% Initialize variables for history
max_steps = 5e6; % Estimate upper bound
timev = zeros(min(max_steps,4096), 1);
history_V = zeros(min(max_steps,4096), 1);
history_x = zeros(min(max_steps,4096), 1);
history_tauf = zeros(min(max_steps,4096), 1);
history_maxgammadot = zeros(min(max_steps,4096), 1);
history_p_max = zeros(min(max_steps,4096), 1); 
history_p_avg = zeros(min(max_steps,4096), 1);
history_p_r = zeros(min(max_steps,4096), 1); 
history_f_r = zeros(min(max_steps,4096), 1);
timev_secs = zeros(min(max_steps,4096), 1); timev_secs(1) = 0;
history_i_export = zeros(min(max_steps,4096), 1); i_export = 0;
%history_p_r = zeros(min(max_steps,4096), 1);

% normalization 
gammadotc = 2 * cs * taub0 / (G * h);

fc = f0 + (a - b) * log(gammadotc/gammadot0);

deltac = rho_c * h / Lambda / fc;

tc = 2 * cs * taub0 / (G * deltac);
%gammac = gammadotc * tc;
gammac = Dc./ h;
if exist('inelast_options','var') && isfield(inelast_options,'gamma_c_override') && ~isempty(inelast_options.gamma_c_override)
    gammac=inelast_options.gamma_c_override;
end
gammac=1; % User-selected manuscript value throughout this suite.

phi = phi_0.*ones(1,ny);
% Clock
tic
save(fullfile(postFolder,'inelast_initialization.mat'),'gammadot_dilatancy_ref','V_slider0','V0','phi_pd','phi_pc','phi_0','p','gammadot','gammac','inelast_reference_profile','-v7');
while current_time < ending_time && count < options.max_steps
    tryagaincount = 0;
    enter = 1;

    % explicit predictor
    V_g = V + (V - V_p) ./ dtp * dt;
    V_mg = 0.5 * (V + V_g);
    gammadot_g = gammadot + (gammadot - gammadot_p) ./ dtp * dt;
    x_g = x + (x - x_p) ./ dtp * dt;
    p_g = p + (p -p_p) ./ dtp * dt;
    %T_g = T + (T - T_p) ./ dtp * dt;

    % previous time step
    V_p = V; % slip rate
    gammadot_p = gammadot; % strain rate
    theta_p = theta; % state
    x_p = x; % localition of the block
    tauf_p = tauf;
    p_p = p;
    T_p = T;
    phi_ep = phi_e;
    phi_pcp = phi_pc;
    phi_pdp = phi_pd;
    phi_pp = phi_p;
    %k_p = k_pnew;
    xi_p = xi;
    current_time_p = current_time;
    count = count + 1;

    while enter || tryagaincount
        switch tryagaincount
            case 0
            % first entry, use the explicitly predicted *_g for computing
            case 1
            % First iteration is considered special and time-step is not
            % reduced, in attempt to obtain a better guess
            p_g = p; T_g = T; x_g = x; V_g = V;
            V_mg = 0.5*(V_g + V_p);
            otherwise
            p_g = 0.5 * (p + p_p); T_g = 0.5 * (T + T_p); x_g = 0.5 * (x + x_p); V_g = 0.5 * (V_p + V);
            V_mg = 0.5 * (V_p + V_g);
        end
        
        enter = 0;
    
        % predictor 
        V_m = 0.5 * (V_g + V_p);
    
        gammadot_m = 0.5 * (gammadot_p + gammadot_g);
    
        theta_predictor = state_theta(theta_p, V_m, dt, Dc);
    
        x_predictor = x_p + V_m * dt;
        
        % corrector 
        [tauf, nb_it] = solvestrength(tauf, sigma_n,  p+pa, theta, k * (Vp(current_time + dt) * (current_time + dt) - x_predictor) + taub0, gammadot0, V0, a, f0, b, Dc, eta, yy_gammadot);
        
        % Reaction kinetics - T in previous step is used in time integration
        decomposition = (1 - xi) .* A .* exp(-Q./ (R.*(T_p + Ta))) .* decomposition_switch; 

        % Dilatancy term
        %n_pl_term = eps_pl *  (1 - V_m .* theta_predictor./ Dc) ./ theta_predictor * 1/beta .* dilatancy;
        %n_pl_term =  eps_pl .* V_m ./ Dc .* (log(theta .* V0 ./ Dc) + log(gammadot ./ gammadot0)) * 1/beta .* dilatancy; % Reference: Eq. 14, Segall et al.,2012
        
        n_pl_term = gammadot_p ./ gammac .* (phi_pdp - eps_pl .* log(max(gammadot_p,realmin)./gammadot_dilatancy_ref)) ./ beta .* dilatancy; % Manuscript rate-only law.
       
        % Implicit solver for p and T in the current step
        %[T, p] = implicit_p_T(T_p, p_p, gammadot, tauf, dt, T_boundarynorth, T_boundarysouth, p_boundarynorth, p_boundarysouth, M_alpha_hy, alpha_th, rho_c, Lambda, dy, n_pl_term, dilatancy, decomposition, rect_m_frac, Pr, Er);

        phi_e = phi_ep +  phi .* (beta_n * (p - p_p) + Lambda_n * (T - T_p));

        phi = phi_e + phi_p;
        
        k_pnew = k_p .* (phi / phi_0).^3;
        

        % Valving 
        M_alpha_hy = [k_pnew(1) 2./(1./k_pnew(1:end-1) + 1./k_pnew(2:end)) k_pnew(end)] ./beta./eta_f .* permeability_evol + M_alpha_hy0 .* (1 - permeability_evol);
        
        %M_alpha_hy = ah_min + (ah_star - ah_min) .* exp(-(sigma_n - p)./sigma_eff_star);
     
       % Update reaction progress xi with T and p at the current step
        if decomposition_switch == 1
          [xi, kDt] = reaction_progress(xi_p, T_p+Ta, T, dt, A, Q, R);
        end

        gammadot = shear_strain(tauf, theta_predictor, gammadot0, V0, f0, b, a, Dc, sigma_n, p+pa);
        
        V = slip_velocity(gammadot, yy_gammadot); 
      
        V_m = 0.5 * (V + V_p); % Mid-point averaging
    
        x = x_p + V_m * dt;  % Slip distance in the current time step
        
        theta = state_theta(theta_p, V_m, dt, Dc);

        phi_pc = phi_pcp - (1 - xi) .* A .* exp(-Q./ (R.*(T + Ta))).* ...
                 beta.* ...
                 decomposition_switch;
        
        phi_pd = (phi_pdp ./ dt + gammadot ./ gammac .* (eps_pl .* log(max(gammadot,realmin)./gammadot_dilatancy_ref))) ./ (1./dt + gammadot./gammac) .* dilatancy;
        %phi_pd = ((V ./Dc .* (phi_0 + eps_pl .* log(gammadot./gammadot0))) + phi_pd./dt) ./ (1./dt + V./Dc) .* ...
        %         dilatancy .* ...
        %        permeability_evol;
        
        phi_p = phi_pd + phi_pc;

        shear_stress = tauf;

        [accuracyabs, ~] = max(abs(p - p_g)./(a.*sigma_n));
    
        if (tryagaincount < miniter || max(abs(p-p_g)./(a.*sigma_n))>0.1*frac || norm(p-p_g)/(norm(p) + 1)>0.1*frac) && tryagaincount<max_iter 
            tryagain = 1;
            tryagaincount = tryagaincount + 1;
            if tryagaincount == 1
                % First iteration refines the prediction and does not halve the
                % time-step. dt = dt;
            else
                dt = dt/2;
            end
        else
            tryagain = 0;
            tryagaincount_hist = tryagaincount;
            tryagaincount = 0;
            if accuracyabs > maxtolvio
                maxtolvio = accuracyabs;
                maxtolviotime = current_time;
            end
        end
        if tryagain == 0 % the time-step has been resolved, update the current total time and next time-step 
            % Check the accepted candidate before real(...) output projections or time advance.
            inelast_candidate=[V(:);tauf(:);p(:);T(:);gammadot(:);theta(:);phi_pd(:);phi_pc(:);phi_p(:);xi(:);M_alpha_hy(:)];
            inelast_audit=inelast_accept_step(inelast_audit,accuracyabs,norm(p-p_g)/(norm(p)+1), ...
                inelast_candidate,V,gammadot,theta,sigma_n-pa-p,tryagaincount_hist,current_time_p+dt,postFolder);
            current_time = current_time_p + dt;
            current_time_year = current_time./secs_per_year;
            
            C1 = 1/16;
            dt_theta1 = C1 * Dc ./ V; % optional: Dc / max(gammadot)
        
            dtp = dt;

            dt = min([dt_theta1 dt_limit], [], 'all');
            if dt > maxdtfac * dtp
                dt = maxdtfac * dtp;
            end
        end
    end
    timev_secs(count+1) = current_time;

    %% Calculate representative values
    % representative state
    % theta_rp = theta_r;
    % theta_r = theta_rp * exp(V * dt / Dc) + Dc / V * (1 - exp(V * dt / Dc));
    % representative pore pressure
    p_r = sigma_n - shear_stress / (f0 + a * log(V_m/V0) + b * log(theta_r / theta0)) - pa;
    % representative f
    f_r = shear_stress / (sigma_n - p_r - pa);
    gamma = gamma + gammadot * dt;
    % save ouputs
    if (mod(count,25) == 0) || (tryagaincount_hist > 2)
    i_export = i_export+1;
    timev(i_export) = real(current_time_year);
    history_V(i_export) = real(V);
    history_x(i_export) = real(x);
    history_tauf(i_export) = real(shear_stress);
    inelast_row.history_state = real(theta);
    inelast_row.history_p = real(p);
    inelast_row.history_gammadot = real(gammadot); gammadot = real(gammadot);
    history_maxgammadot(i_export,:) = max(gammadot);
    inelast_row.history_T = T;
    inelast_row.history_xi = xi;
    history_p_max(i_export) = max(p);
    history_p_avg(i_export) = mean(p);
    inelast_row.history_phi = phi;
    inelast_row.history_phi_e = phi_e;
    inelast_row.history_phi_p = phi_p;
    inelast_row.history_alphahy = M_alpha_hy;
    history_p_r(i_export) = p_r;
    history_f_r(i_export) = f_r;
    history_i_export(i_export) = i_export;
    inelast_stream_fields('append',fullfile(postFolder,'fault_slip_field.mat'),inelast_row,i_export);
    end
    if real(V) > 100
        disp("Model Explosion");
        break;
    end

   
end

% Time cost
t_wallclock = toc;


% Initialize the structure 'writeout' with the required fields
timev               = timev(1:i_export);
history_V           = history_V(1:i_export);
history_x           = history_x(1:i_export);
history_tauf        = history_tauf(1:i_export);
history_maxgammadot = history_maxgammadot(1:i_export);
run_time            = t_wallclock;
history_p_r         = history_p_r(1:i_export);
history_f_r         = history_f_r(1:i_export);
history_p_max         = history_p_max(1:i_export);

% Initialize the structure 'vis' with the required fields

% Save the structure 'writeout' into a MAT-file named 'writeout.mat'
save(fullfile(postFolder,'fault_slip_data.mat'), 'timev', 'history_V', 'history_f_r', 'history_p_max','history_x', 'history_tauf', 'history_maxgammadot', 'run_time', 'a', 'f0' , 'Lambda', 'alpha_th', 'alpha_hy', 'rho_c', 'ny', 'h', 'yy', 'gammadotc', 'deltac', 'fc', 'tc', '-v7');
inelast_stream_report=inelast_stream_fields('finish',fullfile(postFolder,'fault_slip_field.mat'));
assert(inelast_stream_report.rows==i_export);


disp(strcat("Simulation Done after ", num2str(t_wallclock), " seconds."))

% Additional isolated provenance/output; numerical state above is untouched.
inelast_completed=isfinite(current_time) && current_time>=ending_time && ...
    all(isfinite([V(:);p(:);T(:);gammadot(:);phi_pd(:)]));
save(fullfile(postFolder,'inelast_completion.mat'),'inelast_completed','current_time','ending_time','V','p','T','gammadot','phi_pd','phi_pc','phi_p','phi_0','gammadot_dilatancy_ref','gammac','V0','V_slider0','inelast_reference_profile','theta','sigma_n','pa','maxtolvio','inelast_audit','-v7');

save(fullfile(postFolder,'inelast_acceptance_audit.mat'),'inelast_audit','-v7');

result=struct('completed',inelast_completed,'final_days',current_time/86400, ...
    'run_time',run_time,'accepted_steps',count,'samples',numel(timev));
end
