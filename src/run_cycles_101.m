function result=run_cycles_101(options)
inelast_options=options;
%% README
% Current manuscript case cycles_101. Software citation: Yuhan Wang.
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

% Function wrapper: preserve original stepping, refuse existing results.
if nargin < 1, options = struct(); end
if ~isfield(options,'lambda_n'), options.lambda_n=1e-5; end
if ~isfield(options,'stop_after_first_event'), options.stop_after_first_event=false; end
if ~isfield(options,'field_interval'), options.field_interval=1e-4; end
assert(ismember(options.lambda_n,[1e-5,2.4e-5,-2.94e-4]),'Only target or archived-control lambda_n is supported.');
if ~isfield(options,'max_steps'), options.max_steps=Inf; end
if ~isfield(options,'output_directory'), error('An explicit output directory is required'); end
postFolder=char(options.output_directory);
if isfolder(postFolder), error('Output directory already exists: %s',postFolder); end
mkdir(postFolder);
helper_path=fileparts(mfilename('fullpath'));
addpath(helper_path);
cleanup_path=onCleanup(@() rmpath(helper_path));

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
permeability_evol = 1;
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
% Sole active parameter change: retain Lambda, storage, and all other inputs.
% Lambda_f is implied by these independent inputs and is not read in stepping.
if options.lambda_n > 0
    Lambda_n = options.lambda_n;
    Lambda_f = Lambda_n + Lambda * beta / phi_0;
end


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

% Bounded diagnostic storage; none of these quantities enter stepping.
history=zeros(65536,20); i_export=0;
peak_V=0; peak_time=0; peak_state=struct();
minimum_sigma_e=Inf; maximum_dilatancy_step=0;
maximum_pressure_retry=0; maximum_force_residual=0;
minimum_absolute_pressure=Inf; minimum_porosity=Inf; maximum_porosity=-Inf;
minimum_alpha=Inf; maximum_alpha=-Inf; maximum_p=-Inf; minimum_p=Inf; maximum_T=-Inf;
maximum_boundary_p=0; maximum_boundary_T=0;
event_count=0; event_active=false; last_progress_wall=0;
field_count=0; field_active=false; field_complete=false; field_last_time=-Inf;
field_interval=options.field_interval; field_capacity=2048; field_max=50000;
spatial=struct();
spatial_names={'p','T','gammadot','xi','phi_e','phi_p','phi_p_updated','phi_pd','phi_pc','phi'};
for j=1:numel(spatial_names),spatial.(spatial_names{j})=zeros(field_capacity,ny,'single');end
spatial.alphahy=zeros(field_capacity,ny+1,'single');
scalar_field_names={'time_seconds','V','tauf','x','theta','force_residual'};
for j=1:numel(scalar_field_names),spatial.(scalar_field_names{j})=zeros(field_capacity,1);end
field_stop_time=Inf;

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
solver_clock=tic;
parameters=struct('case_code','101','gamma_c',gammac,'D_RS',Dc,'h',h,'ny',ny, ...
 'a',a,'b',b,'f0',f0,'V0',V0,'loading_rate',Vp(1),'k',k,'eta',eta, ...
 'sigma_n',sigma_n,'pa',pa,'Ta',Ta,'Lambda',Lambda,'rho_c',rho_c, ...
 'alpha_th',alpha_th,'alpha_hy',alpha_hy,'eps_pl',eps_pl,'beta',beta, ...
 'dilatancy',dilatancy,'permeability_evol',permeability_evol, ...
 'decomposition_switch',decomposition_switch,'rect_m_frac',rect_m_frac, ...
 'random_seed',334500,'ending_time',ending_time,'dt_limit',dt_limit, ...
 'frac',frac,'maxdtfac',maxdtfac,'state_step_fraction',1/16, ...
 'boundary_condition',options.boundary_type,'state_law','slip aging','friction_form','archive sinh', ...
 'phi0',phi_0,'Sc',beta,'beta_n',beta_n,'beta_f',beta_f,'lambda_n',Lambda_n, ...
 'lambda_f',Lambda_f,'eta_f',eta_f,'k0',k_p,'Pr',Pr,'Er',Er,'rho_f',rhof, ...
 'A',A,'Q',Q,'R',R,'mbar',rect_m_frac,'thermal_temperature_offset',Ta, ...
 'porosity_model','unaltered archived incremental and effective inelastic update', ...
 'temperature_convention','unaltered archived Ta=170 and reaction_progress(Tp+Ta,T)', ...
 'field_phi_e_convention','elastic change = legacy phi_e minus phi0', ...
 'field_phi_p_convention','lagged effective phi_p actually used in permeability; phi_p_updated is subsequent update');
initial_fields=struct('p',p,'T',T,'xi',xi,'phi_e',zeros(1,ny),'phi_p',zeros(1,ny), ...
 'phi_e_legacy',phi_e,'phi_p_updated',phi_p,'phi0',phi_0);
save(fullfile(postFolder,'parameters.mat'),'parameters','initial_fields');
i_export=1;
history(1,:)=[0,V,x,tauf,max(gammadot),max(p),f0,0,max(T),max(phi_pd), ...
 mean(xi),max(xi),min(p),phi_0,phi_0,1,1,k*(Vp(0)*0-x)+taub0-tauf-eta*V, ...
 trapz(yy,xi)/h,theta(1)];
fprintf('START case=%s gamma_c=%g horizon_days=%g\n',parameters.case_code,gammac,ending_time/86400);
try
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
        [T, p] = implicit_p_T(T_p, p_p, gammadot, tauf, dt, T_boundarynorth, T_boundarysouth, p_boundarynorth, p_boundarysouth, M_alpha_hy, alpha_th, rho_c, Lambda, dy, n_pl_term, dilatancy, decomposition, rect_m_frac, Pr, Er, release_boundary);

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
    % Passive checks never modify a state or relax an archive tolerance.
    if ~isreal([V,tauf,p,T,gammadot,theta,phi_pd,phi_pc,phi_e,phi,xi,M_alpha_hy]) || ...
       any(~isfinite([V,tauf,p,T,gammadot,theta,phi_pd,phi_pc,phi_e,phi,xi,M_alpha_hy])) || ...
       min(sigma_n-pa-p)<=0 || V<=0 || any(gammadot<=0) || any(theta<=0) || ...
       any(phi<=0) || any(M_alpha_hy<=0)
        error('rerun:InvalidState','Nonphysical accepted state at %.16g seconds',current_time);
    end
    if tryagaincount_hist>=max_iter && (accuracyabs>0.1*frac || norm(p-p_g)/(norm(p)+1)>0.1*frac)
        error('rerun:PressureRetries','Pressure acceptance failed at retry limit');
    end
    minimum_sigma_e=min(minimum_sigma_e,min(sigma_n-pa-p));
    minimum_absolute_pressure=min(minimum_absolute_pressure,min(pa+p));
    minimum_porosity=min(minimum_porosity,min(phi));maximum_porosity=max(maximum_porosity,max(phi));
    minimum_alpha=min(minimum_alpha,min(M_alpha_hy));maximum_alpha=max(maximum_alpha,max(M_alpha_hy));
    maximum_p=max(maximum_p,max(p));minimum_p=min(minimum_p,min(p));maximum_T=max(maximum_T,max(T));
    maximum_boundary_p=max(maximum_boundary_p,max(abs(p([1,end]))));
    maximum_boundary_T=max(maximum_boundary_T,max(abs(T([1,end]))));
    maximum_dilatancy_step=max(maximum_dilatancy_step,dtp*max(gammadot)/gammac);
    maximum_pressure_retry=max(maximum_pressure_retry,tryagaincount_hist);
    force_residual=k*(Vp(current_time)*current_time-x)+taub0-tauf-eta*V;
    maximum_force_residual=max(maximum_force_residual,abs(force_residual));
    if V>peak_V
        peak_V=V; peak_time=current_time;
        peak_state=struct('time',current_time,'V',V,'p',p,'T',T,'gammadot',gammadot, ...
          'theta',theta,'tauf',tauf,'x',x,'phi',phi,'phi_e_legacy',phi_e, ...
          'phi_p_used',phi-phi_e,'phi_p_updated',phi_p,'phi_pc',phi_pc,'phi_pd',phi_pd, ...
          'xi',xi,'alphahy',M_alpha_hy,'force_residual',force_residual);
    end
    if ~event_active && V>=0.1, event_active=true; end
    if event_active && V<1e-4, event_count=event_count+1; event_active=false; end
    if ~field_active && ~field_complete && V>=1e-4, field_active=true; end
    if field_active && ~field_complete && current_time-field_last_time>=field_interval
        if field_count==field_max
            keep=1:2:field_count; nkeep=numel(keep); names=fieldnames(spatial);
            for j=1:numel(names),spatial.(names{j})(1:nkeep,:)=spatial.(names{j})(keep,:);end
            field_count=nkeep;field_interval=2*field_interval;
        elseif field_count==field_capacity
            field_capacity=min(field_max,2*field_capacity);names=fieldnames(spatial);
            for j=1:numel(names),spatial.(names{j})(field_capacity,end)=0;end
        end
        field_count=field_count+1;
        spatial.time_seconds(field_count)=current_time;spatial.V(field_count)=V;
        spatial.tauf(field_count)=tauf;spatial.x(field_count)=x;spatial.theta(field_count)=theta(1);
        spatial.force_residual(field_count)=force_residual;
        spatial.p(field_count,:)=single(p);spatial.T(field_count,:)=single(T);
        spatial.gammadot(field_count,:)=single(gammadot);spatial.alphahy(field_count,:)=single(M_alpha_hy);
        spatial.xi(field_count,:)=single(xi);spatial.phi_e(field_count,:)=single(phi_e-phi_0);
        spatial.phi_p(field_count,:)=single(phi-phi_e);spatial.phi_p_updated(field_count,:)=single(phi_p);
        spatial.phi_pd(field_count,:)=single(phi_pd);spatial.phi_pc(field_count,:)=single(phi_pc);
        spatial.phi(field_count,:)=single(phi);
        field_last_time=current_time;
        if event_count>=1 && isinf(field_stop_time),field_stop_time=current_time+10;end
        if current_time>=field_stop_time,field_complete=true;end
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
    % Dense diagnostic sampling is separate from the original state projection.
    legacy_projection_due=(mod(count,25)==0) || (tryagaincount_hist>2);
    if legacy_projection_due || V>=1e-4 || current_time>=ending_time
        if i_export==size(history,1),history(2*size(history,1),20)=0;end
        i_export=i_export+1;
        history(i_export,:)=[real(current_time_year),real(V),real(x),real(tauf), ...
          max(real(gammadot)),max(p),f_r,dtp,max(T),max(phi_pd),mean(xi),max(xi), ...
          min(p),min(phi),max(phi),min(M_alpha_hy)/alpha_hy,max(M_alpha_hy)/alpha_hy, ...
          force_residual,trapz(yy,xi)/h,theta(1)];
    end
    if legacy_projection_due,gammadot=real(gammadot);end
    elapsed=toc(solver_clock);
    if elapsed-last_progress_wall>=30
        fprintf('PROGRESS case=%s steps=%d days=%.9g V=%.8g events=%d wall_s=%.2f samples=%d\n', ...
            parameters.case_code,count,current_time/86400,V,event_count,elapsed,i_export);
        save(fullfile(postFolder,'progress.mat'),'parameters','count','current_time','V','event_count','peak_V','elapsed');
        last_progress_wall=elapsed;
    end
    if real(V)>100,error('rerun:Explosion','Slip rate exceeds archive 100 m/s guard');end
    if options.stop_after_first_event && field_complete,break,end
end
catch failure
    save(fullfile(postFolder,'failure.mat'),'failure','parameters','count','current_time','dt','dtp', ...
      'V','x','tauf','theta','p','T','gammadot','xi','phi_pd','phi_p','phi_e','M_alpha_hy','history','i_export');
    rethrow(failure);
end

run_time=toc(solver_clock);
history=history(1:i_export,:);
timev=history(:,1); history_V=history(:,2); history_x=history(:,3); history_tauf=history(:,4);
history_maxgammadot=history(:,5); history_p_max=history(:,6); history_f_r=history(:,7);
history_dt=history(:,8); history_T_max=history(:,9); history_phi_pd_max=history(:,10);
history_xi_mean=history(:,11); history_xi_max=history(:,12);
completed=current_time>=ending_time;accepted_steps=count;final_time=current_time;
pilot_complete=options.stop_after_first_event && field_complete;
history_p_min=history(:,13);history_phi_min=history(:,14);history_phi_max=history(:,15);
history_alpha_min=history(:,16);history_alpha_max=history(:,17);history_force_residual=history(:,18);
history_xi_mean_weighted=history(:,19);history_theta=history(:,20);
gamma_c=gammac;
save(fullfile(postFolder,'inelast_acceptance_audit.mat'),'inelast_audit','-v7');
diagnostics=struct('minimum_sigma_e',minimum_sigma_e,'maximum_dilatancy_step',maximum_dilatancy_step, ...
 'maximum_pressure_retry',maximum_pressure_retry,'maximum_force_residual_Pa',maximum_force_residual, ...
 'peak_V',peak_V,'peak_time',peak_time,'event_count',event_count,'event_active',event_active, ...
 'fields_complete',field_complete,'maximum_boundary_p',maximum_boundary_p, ...
 'maximum_boundary_T',maximum_boundary_T,'minimum_absolute_pressure',minimum_absolute_pressure, ...
 'minimum_phi',minimum_porosity,'maximum_phi',maximum_porosity, ...
 'minimum_alpha_ratio',minimum_alpha/alpha_hy,'maximum_alpha_ratio',maximum_alpha/alpha_hy, ...
 'maximum_p',maximum_p,'minimum_p',minimum_p,'maximum_T',maximum_T);
save(fullfile(postFolder,'fault_slip_data.mat'),'timev','history_V','history_f_r','history_p_max','history_x', ...
 'history_tauf','history_maxgammadot','run_time','a','b','V0','f0','Lambda','alpha_th','alpha_hy','rho_c', ...
 'ny','h','yy','gammadotc','deltac','fc','tc','gamma_c','parameters','completed','accepted_steps','final_time', ...
 'history_dt','history_T_max','history_phi_pd_max','history_xi_mean','history_xi_max','diagnostics', ...
 'history_p_min','history_phi_min','history_phi_max','history_alpha_min','history_alpha_max', ...
 'history_force_residual','history_xi_mean_weighted','history_theta','pilot_complete','-v7');
names=fieldnames(spatial);
for j=1:numel(names),spatial.(names{j})=spatial.(names{j})(1:field_count,:);end
spatial.yy=yy;spatial.h=h;spatial.interval=field_interval;spatial.complete=field_complete;
spatial.phi_e_initial=zeros(1,ny);spatial.phi_p_initial=zeros(1,ny);
spatial.phi_e_convention='signed elastic change relative to initial phi0';
spatial.phi_p_convention='actual lagged effective inelastic porosity in permeability';
save(fullfile(postFolder,'first_event_fields.mat'),'spatial','parameters','-v7');
save(fullfile(postFolder,'final_state.mat'),'parameters','completed','count','current_time','dt','dtp', ...
 'V','x','tauf','theta','p','T','gammadot','xi','phi_pd','phi_pc','phi_p','phi_e','phi', ...
 'M_alpha_hy','gamma','peak_state','initial_fields','pilot_complete');
result=struct('completed',completed,'pilot_complete',pilot_complete,'final_days',current_time/86400,'run_time',run_time, ...
 'accepted_steps',count,'samples',i_export,'diagnostics',diagnostics);
fprintf('DONE case=%s completed=%d days=%.12g steps=%d wall_s=%.3f peak_V=%.9g events=%d\n', ...
 parameters.case_code,completed,current_time/86400,count,run_time,peak_V,event_count);
end
