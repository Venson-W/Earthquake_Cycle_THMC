function result=run_dense_001(options)
inelast_options=options;
%% README
% Current manuscript case dense_001. Software citation: Yuhan Wang.
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

% Safe wrapper: isolated output, exact archived material parameters and laws.
if nargin<1, options=struct(); end
if ~isfield(options,'output_directory'),error('An explicit fresh output directory is required.');end
if ~isfield(options,'event_max_dt'),options.event_max_dt=2e-5;end
% Capture follows this run's accepted trajectory; no historical timing input.
if isfield(options,'capture_start_seconds') && isfinite(options.capture_start_seconds)
    error('dense:ExternalCaptureTime','Dense capture must be derived from its own accepted state.');
end
options.capture_start_seconds=NaN;
options.capture_fraction=0.1;
options.capture_mode='accepted_state_threshold';
if ~isfield(options,'max_fields'),options.max_fields=400000;end
if ~isfield(options,'max_steps'),options.max_steps=Inf;end
assert(isscalar(options.event_max_dt) && options.event_max_dt>0,'event_max_dt must be positive.');
postFolder=char(options.output_directory);
if isfolder(postFolder) || isfile(postFolder),error('Output exists: %s',postFolder);end
mkdir(postFolder);
helper_path=fileparts(mfilename('fullpath'));
addpath(helper_path);
cleanup_path=onCleanup(@()rmpath(helper_path));

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
dilatancy = 1;

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

% normalization 
gammadotc = 2 * cs * taub0 / (G * h);

fc = f0 + (a - b) * log(gammadotc/gammadot0);

deltac = rho_c * h / Lambda / fc;

tc = 2 * cs * taub0 / (G * deltac);
%gammac = gammadotc * tc;
gammac = Dc./ (h);
if exist('inelast_options','var') && isfield(inelast_options,'gamma_c_override') && ~isempty(inelast_options.gamma_c_override)
    gammac=inelast_options.gamma_c_override;
end
gammac=1; % User-selected manuscript value throughout this suite.

phi = phi_0.*ones(1,ny);
% Passive storage is independent of the original projection/export cadence.
history=zeros(65536,19);i_export=0;
field_capacity=8192;field_count=0;
spatial=struct('time_seconds',zeros(field_capacity,1),'V',zeros(field_capacity,1), ...
 'step_seconds',zeros(field_capacity,1),'p',zeros(field_capacity,ny,'single'), ...
 'T',zeros(field_capacity,ny,'single'),'gammadot',zeros(field_capacity,ny,'single'), ...
 'xi',zeros(field_capacity,ny,'single'));
localization_threshold=pi^2*a/(f0*(f0-2*a))*rho_c/Lambda*(alpha_hy+alpha_th)/h;
capture_threshold=options.capture_fraction*localization_threshold;capture_started=false;
capture_trigger_step_seconds=NaN;maximum_dense_sample_interval=0;
assert(capture_threshold>V && capture_threshold<localization_threshold,'Capture trigger must precede localization.');
first_localization_time=NaN;rapid_event_seen=false;rapid_event_end=NaN;
field_complete=false;last_progress_wall=0;
maximum_force_residual=0;maximum_pressure_retry=0;
minimum_sigma_e=Inf;minimum_absolute_pressure=Inf;minimum_step=Inf;maximum_capture_step=0;
peak_V=0;peak_time=NaN;peak_p=-Inf;peak_T=-Inf;
parameters=struct('case_code','001','h',h,'ny',ny,'a',a,'b',b,'D_RS',Dc, ...
 'f0',f0,'V0',V0,'loading_rate',Vp(1),'k',k,'eta',eta, ...
 'sigma_n',sigma_n,'pa',pa,'Ta',Ta,'Lambda',Lambda,'rho_c',rho_c, ...
 'alpha_th',alpha_th,'alpha_hy',alpha_hy,'eps_pl',eps_pl,'Sc',beta,'beta',beta, ...
 'dilatancy',dilatancy,'permeability_evol',permeability_evol, ...
 'decomposition_switch',decomposition_switch,'rect_m_frac',rect_m_frac,'mbar',rect_m_frac, ...
 'random_seed',334500,'ending_time',ending_time,'dt_limit',dt_limit,'frac',frac, ...
 'maxdtfac',maxdtfac,'state_step_fraction',1/16,'boundary_condition',options.boundary_type, ...
 'phi0',phi_0,'beta_n',beta_n,'beta_f',beta_f,'lambda_n',Lambda_n,'lambda_f',Lambda_f, ...
 'eta_f',eta_f,'k0',k_p,'Pr',Pr,'Er',Er,'A',A,'Q',Q,'R',R, ...
 'event_max_dt',options.event_max_dt,'capture_start_seconds',options.capture_start_seconds, ...
 'capture_mode',options.capture_mode,'capture_fraction',options.capture_fraction,'capture_threshold',capture_threshold, ...
 'localization_threshold',localization_threshold,'field_sampling','every accepted step; no thinning', ...
 'first_event_stop','V reached 0.1 then below 1e-4 plus 0.05s, and 0.30s after first localization crossing', ...
 'temperature_convention','unaltered archive: Ta=170 and reaction_progress(Tp+Ta,T)', ...
 'porosity_model','exact original case-specific dormant diagnostic law');
initial_fields=struct('p',p,'T',T,'xi',xi,'gammadot',gammadot,'theta',theta, ...
 'V',V,'tauf',tauf,'x',x,'phi_e',phi_e,'phi_p',phi_p,'phi_pd',phi_pd,'phi_pc',phi_pc);
save(fullfile(postFolder,'parameters.mat'),'parameters','initial_fields','options');
solver_clock=tic;
fprintf('START case=%s event_max_dt=%.8g capture_start=%.12g\n',parameters.case_code,options.event_max_dt,options.capture_start_seconds);
try
save(fullfile(postFolder,'inelast_initialization.mat'),'gammadot_dilatancy_ref','V_slider0','V0','phi_pd','phi_pc','phi_0','p','gammadot','gammac','inelast_reference_profile','-v7');

while current_time < ending_time && count < options.max_steps
    % Preserve archive steps up to the accepted trigger. Refine only later steps.
    if capture_started && isfinite(options.event_max_dt)
        dt=min(dt,options.event_max_dt);
    end

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
        [T, p] = implicit_p_T(T_p, p_p, gammadot, tauf, dt, T_boundarynorth, T_boundarysouth, p_boundarynorth, p_boundarysouth, M_alpha_hy, alpha_th, rho_c, Lambda, dy, n_pl_term, dilatancy, decomposition, rect_m_frac, Pr, Er, release_boundary);

        phi_e = phi_0 .* (beta_n * p + Lambda_n * T);

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

    %% Calculate representative values
    % representative state
    % theta_rp = theta_r;
    % theta_r = theta_rp * exp(V * dt / Dc) + Dc / V * (1 - exp(V * dt / Dc));
    % representative pore pressure
    p_r = sigma_n - shear_stress / (f0 + a * log(V_m/V0) + b * log(theta_r / theta0)) - pa;
    % representative f
    f_r = shear_stress / (sigma_n - p_r - pa);
    gamma = gamma + gammadot * dt;
    % Original scalar projection remains at its original accepted-step cadence.
    legacy_export_due=(mod(count,25)==0) || (tryagaincount_hist>2);
    if legacy_export_due,gammadot=real(gammadot);end
    accepted_state=[V,tauf,p,T,gammadot,theta,xi,M_alpha_hy];
    if ~isreal(accepted_state) || any(~isfinite(accepted_state)) || ...
       V<=0 || any(gammadot<=0) || any(theta<=0) || min(sigma_n-pa-p)<=0
        error('dense:InvalidState','Invalid accepted state at %.16g seconds',current_time);
    end
    if tryagaincount_hist>=max_iter && (accuracyabs>0.1*frac || norm(p-p_g)/(norm(p)+1)>0.1*frac)
        error('dense:PressureRetries','Pressure acceptance reached failed retry limit');
    end
    force_residual=k*(Vp(current_time)*current_time-x)+taub0-tauf-eta*V;
    maximum_force_residual=max(maximum_force_residual,abs(force_residual));
    maximum_pressure_retry=max(maximum_pressure_retry,tryagaincount_hist);
    minimum_sigma_e=min(minimum_sigma_e,min(sigma_n-pa-p));
    minimum_absolute_pressure=min(minimum_absolute_pressure,min(pa+p));
    minimum_step=min(minimum_step,dtp);
    peak_p=max(peak_p,max(p));peak_T=max(peak_T,max(T));
    if V>peak_V,peak_V=V;peak_time=current_time;end
    if ~capture_started && V>=capture_threshold
        if V>=localization_threshold
            error('dense:MissedCapture','First accepted capture trigger already reached localization.');
        end
        capture_started=true;options.capture_start_seconds=current_time;
        capture_trigger_step_seconds=dtp;
        parameters.capture_start_seconds=current_time;
        parameters.capture_trigger_slip_rate=V;
        parameters.capture_trigger_step_seconds=dtp;
        save(fullfile(postFolder,'parameters.mat'),'parameters','initial_fields','options');
    end
    if isnan(first_localization_time) && V>=localization_threshold,first_localization_time=current_time;end
    if V>=0.1,rapid_event_seen=true;end
    if rapid_event_seen && V<1e-4 && isnan(rapid_event_end),rapid_event_end=current_time;end
    capture_due=capture_started;
    if legacy_export_due || capture_due
        if i_export==size(history,1),history(2*size(history,1),19)=0;end
        i_export=i_export+1;
        history(i_export,:)=[real(current_time_year),real(V),real(x),real(shear_stress), ...
          max(real(gammadot)),max(p),real(f_r),dtp,max(T),mean(xi),max(xi),min(p), ...
          force_residual,theta(1),count,legacy_export_due,tryagaincount_hist,real(p_r),mean(p)];
    end
    if capture_due
        if field_count>=options.max_fields,error('dense:FieldCapacity','Field limit reached without thinning.');end
        if field_count==field_capacity
            field_capacity=min(options.max_fields,2*field_capacity);
            names=fieldnames(spatial);
            for j=1:numel(names),spatial.(names{j})(field_capacity,end)=0;end
        end
        if field_count>0
            dense_sample_interval=current_time-spatial.time_seconds(field_count);
            maximum_dense_sample_interval=max(maximum_dense_sample_interval,dense_sample_interval);
            if dense_sample_interval>options.event_max_dt+16*eps(max(1,abs(current_time)))
                error('dense:CaptureStep','A dense sample interval exceeds the requested event timestep.');
            end
        end
        field_count=field_count+1;
        spatial.time_seconds(field_count)=current_time;spatial.V(field_count)=V;
        spatial.step_seconds(field_count)=dtp;
        spatial.p(field_count,:)=single(p);spatial.T(field_count,:)=single(T);
        spatial.gammadot(field_count,:)=single(gammadot);spatial.xi(field_count,:)=single(xi);
        if field_count>1,maximum_capture_step=max(maximum_capture_step,dtp);end
    end
    if ~isnan(rapid_event_end) && ~isnan(first_localization_time) && ...
       current_time>=rapid_event_end+0.05 && current_time>=first_localization_time+0.30
        field_complete=true;
    end
    elapsed=toc(solver_clock);
    if elapsed-last_progress_wall>=30
        fprintf('PROGRESS case=%s steps=%d seconds=%.12g V=%.7g fields=%d wall_s=%.2f\n', ...
          parameters.case_code,count,current_time,V,field_count,elapsed);
        save(fullfile(postFolder,'progress.mat'),'parameters','count','current_time','V','field_count','elapsed');
        last_progress_wall=elapsed;
    end
    if real(V)>100,error('dense:Explosion','Slip rate exceeds archived 100m/s guard');end
    if field_complete,break,end
end

catch failure
    save(fullfile(postFolder,'failure.mat'),'failure','parameters','count','current_time','dt','dtp', ...
      'V','x','tauf','theta','p','T','gammadot','xi','phi_pd','phi_pc','phi_p','phi_e', ...
      'M_alpha_hy','i_export','field_count');
    rethrow(failure);
end
run_time=toc(solver_clock);
history=history(1:i_export,:);
timev=history(:,1);history_V=history(:,2);history_x=history(:,3);history_tauf=history(:,4);
history_maxgammadot=history(:,5);history_p_max=history(:,6);history_f_r=history(:,7);
history_dt=history(:,8);history_T_max=history(:,9);history_xi_mean=history(:,10);history_xi_max=history(:,11);
history_p_min=history(:,12);history_force_residual=history(:,13);history_theta=history(:,14);
history_step_count=history(:,15);legacy_export_flag=logical(history(:,16));
history_retry_count=history(:,17);history_p_r=history(:,18);history_p_avg=history(:,19);
accepted_steps=count;final_time=current_time;completed=current_time>=ending_time;
first_event_complete=field_complete;
save(fullfile(postFolder,'inelast_acceptance_audit.mat'),'inelast_audit','-v7');
diagnostics=struct('maximum_force_residual_Pa',maximum_force_residual, ...
 'maximum_pressure_retry',maximum_pressure_retry,'minimum_sigma_e',minimum_sigma_e, ...
 'minimum_absolute_pressure',minimum_absolute_pressure,'minimum_step_seconds',minimum_step, ...
 'maximum_capture_step_seconds',maximum_capture_step,'peak_V',peak_V,'peak_time',peak_time, ...
 'capture_start_seconds',options.capture_start_seconds,'capture_trigger_threshold',capture_threshold, ...
 'capture_trigger_step_seconds',capture_trigger_step_seconds,'maximum_dense_sample_interval',maximum_dense_sample_interval, ...
 'maximum_p',peak_p,'maximum_T',peak_T,'first_localization_time',first_localization_time, ...
 'rapid_event_end',rapid_event_end,'field_count',field_count,'first_event_complete',first_event_complete);
save(fullfile(postFolder,'fault_slip_data.mat'),'timev','history_V','history_f_r','history_p_max', ...
 'history_x','history_tauf','history_maxgammadot','run_time','a','b','V0','f0','Lambda', ...
 'alpha_th','alpha_hy','rho_c','ny','h','yy','gammadotc','deltac','fc','tc','parameters', ...
 'completed','first_event_complete','accepted_steps','final_time','diagnostics','history_dt', ...
 'history_T_max','history_xi_mean','history_xi_max','history_p_min','history_force_residual', ...
 'history_theta','history_step_count','legacy_export_flag','history_retry_count', ...
 'history_p_r','history_p_avg','-v7');
names=fieldnames(spatial);
for j=1:numel(names),spatial.(names{j})=spatial.(names{j})(1:field_count,:);end
spatial.yy=yy;spatial.h=h;spatial.complete=field_complete;spatial.parameters=parameters;
spatial.field_sampling='Every accepted state from own-trajectory trigger; no thinning. Trigger sample ends an inherited coarse step; subsequent sample intervals use event_max_dt.';
spatial.first_localization_time=first_localization_time;spatial.rapid_event_end=rapid_event_end;
spatial.diagnostics=diagnostics;
save(fullfile(postFolder,'first_event_fields.mat'),'spatial','parameters','-v7.3');
% Full current and previous dynamic state: sufficient to inspect the final update.
save(fullfile(postFolder,'final_state.mat'),'parameters','options','initial_fields','completed', ...
 'first_event_complete','count','current_time','current_time_p','dt','dtp','V','V_p','V_g','V_m', ...
 'x','x_p','x_g','x_predictor','tauf','tauf_p','theta','theta_p','theta_predictor','theta_r', ...
 'p','p_p','p_g','T','T_p','T_g','gammadot','gammadot_p','gammadot_g','gammadot_m', ...
 'xi','xi_p','phi_pd','phi_pdp','phi_pc','phi_pcp','phi_p','phi_e','phi_ep','phi', ...
 'M_alpha_hy','M_alpha_hy0','k_pnew','gamma','n_pl_term','decomposition', ...
 'shear_stress','p_r','f_r','maxtolvio','diagnostics');
result=struct('first_event_complete',field_complete,'final_seconds',current_time, ...
 'run_time',run_time,'accepted_steps',count,'fields',field_count,'diagnostics',diagnostics);
fprintf('DONE case=%s first_event_complete=%d seconds=%.12g steps=%d fields=%d wall_s=%.3f\n', ...
 parameters.case_code,field_complete,current_time,count,field_count,run_time);
end
