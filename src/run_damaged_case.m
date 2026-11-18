function result = run_damaged_case(cfg)
%RUN_DAMAGED_CASE Current MRSF solver with geometric chemical porosity.
% Yuhan Wang, ORCID 0000-0003-0490-5977. See CITATION.cff.
% All outputs are new files in cfg.output_directory; archives are untouched.
% p and T are changes from pa and T0; phi_e is an elastic porosity CHANGE.
if isfield(cfg,'selftest') && cfg.selftest
 result=selftest(cfg);return
end
cfg = defaults(cfg);
cfg.boundary_type=validatestring(cfg.boundary_type,{'dirichlet','robin'});
out = cfg.output_directory;
assert(~isfolder(out) && ~isfile(out), 'Output already exists; choose a new directory.');
mkdir(out);
rng(cfg.random_seed);
n=cfg.ny; h=cfg.h; yy=linspace(-h/2,h/2,n); dy=h/(n-1);
weights=dy*ones(1,n); weights([1,end])=dy/2;
cfg.weights=weights; cfg.yy=yy;
cfg.permeability_evol=~strcmp(cfg.case_code,'0') && cfg.case_code(1)=='1';
cfg.decomposition_switch=~strcmp(cfg.case_code,'0') && cfg.case_code(2)=='1';
cfg.dilatancy=~strcmp(cfg.case_code,'0') && cfg.case_code(3)=='1';
cfg.thermal_pressurization=~strcmp(cfg.case_code,'0');
cfg.gammadot0=cfg.V0/h;
% Uniform initial reference precedes intentional pressure perturbations.
if strcmpi(cfg.model,'gaussian')
 cfg.gammadot_dilatancy_ref=cfg.V0*cfg.gaussian_shape;
else
 cfg.gammadot_dilatancy_ref=cfg.gammadot0*ones(1,cfg.ny);
end
assert(cfg.gamma_c==1,'Manuscript gamma_c must equal1.');
cfg.taub0=cfg.f0*(cfg.sigma_n-cfg.pa);
cfg.k=cfg.R_k*(cfg.sigma_n-cfg.pa)*(cfg.b-cfg.a)/cfg.D_RS;
if isfield(cfg,'k_s'),cfg.k=cfg.k_s;end
cfg.robin=thmc_boundary_options(cfg,cfg.alpha_hy,cfg.alpha_th,cfg.rho_c,h);
fixed=struct('ny',n,'alpha_th',cfg.alpha_th,'rho_c',cfg.rho_c, ...
 'Lambda',cfg.Lambda*cfg.thermal_pressurization,'dy',dy,'rect_m_frac',cfg.mbar, ...
 'Pr',cfg.Pr,'Er',cfg.Er,'robin',cfg.robin,'temperature_ambient',0,'pressure_ambient',0);
step_p_T=make_p_T_stepper(fixed);
master_y=linspace(-h/2,h/2,300);master_p=cfg.pa*.001*randn(1,300);
initial_p=interp1(master_y,master_p,yy,'linear');
if strcmpi(cfg.boundary_type,'dirichlet'),initial_p([1,end])=0;end
s=struct('time',0,'p',initial_p,'T',zeros(1,n), ...
 'xi',zeros(1,n),'phi_pd',zeros(1,n),'theta',cfg.D_RS/cfg.V0,'x',0, ...
 'V',cfg.V0,'tauf',cfg.taub0,'gammadot',cfg.gammadot_dilatancy_ref, ...
 'gamma_accumulated',zeros(1,n));
s=porosity(s,cfg);
[s.tauf,s.V,s.gammadot,s.theta,s.x,force]=mechanics(s.p,s,0,cfg,s.tauf);
parameters=cfg; initial_state=s;
save(fullfile(out,'parameters.mat'),'parameters','initial_state');
history=zeros(65536,18); ih=1;
history(ih,:)=[0,s.V,s.x,s.tauf,max(s.gammadot),max(s.p),0,max(s.T), ...
 max(s.phi_pd),sum(weights.*s.xi)/h,max(s.xi),min(s.p),min(s.phi),max(s.phi), ...
 min(s.alphahy)/cfg.alpha_hy,max(s.alphahy)/cfg.alpha_hy,force,s.theta];
fields=struct(); nf=0; capacity=2048; field_interval=cfg.field_interval;
field_names={'p','T','gammadot','alphahy','xi','phi_e','phi_pd','phi_pc','phi_p'};
for j=1:numel(field_names)
 w=n; if strcmp(field_names{j},'alphahy'),w=n+1;end
 fields.(field_names{j})=zeros(capacity,w,'single');
end
fields.time_seconds=zeros(capacity,1); fields.V=zeros(capacity,1); fields.step_seconds=zeros(capacity,1);
fields.tauf=zeros(capacity,1); fields.x=zeros(capacity,1); fields.theta=zeros(capacity,1);
dense_event=false; event_active=false; event_count=0; fields_active=false; fields_complete=false;
field_stop=Inf; last_field=-Inf; dt=1e-7*cfg.dt_scale;
count=0; rejected=0; max_force=abs(force); minimum_sigma_e=Inf;
maximum_coupling_iterations=0; max_reaction_step=0; maximum_phi=0; minimum_phi=Inf;
peak_state=s; peak_V=s.V; maximum_boundary_p=0; maximum_boundary_T=0;
last_progress=0; timer=tic;
fprintf('START case=%s model=%s days=%g ny=%d dt_scale=%g\n',cfg.case_code,cfg.model,cfg.end_seconds/86400,n,cfg.dt_scale);
try
while s.time<cfg.end_seconds && count<cfg.max_steps
 if ~fields_complete && s.V>=cfg.dense_velocity_threshold,dense_event=true;end
 dt=min([dt,cfg.end_seconds-s.time,cfg.dt_limit*cfg.dt_scale]);
 if dense_event && ~fields_complete,dt=min(dt,cfg.first_event_dt_max);end
 accepted=false;
 for retry=0:55
  [candidate,ok,it,force]=coupled_step(s,dt,cfg,step_p_T);
  if ok
   change=max([max(abs(candidate.p-s.p))/(cfg.pressure_step*cfg.dt_scale), ...
       max(abs(candidate.T-s.T))/(cfg.temperature_step*cfg.dt_scale), ...
       max(candidate.xi-s.xi)/(cfg.reaction_step*cfg.dt_scale), ...
       abs(log(max(candidate.V,realmin)/max(s.V,realmin)))/(cfg.log_velocity_step*cfg.dt_scale)]);
   if change<=1.01 || dt<1e-9
    accepted=true; break
   end
  end
  rejected=rejected+1; dt=dt/2;
  if s.time+dt==s.time, error('damaged:TimeStepUnderflow','Time step underflow at %.16g s.',s.time); end
 end
 if ~accepted, error('damaged:StepFailure','Failed coupled step at %.16g s.',s.time); end
 candidate.gamma_accumulated=s.gamma_accumulated+.5*dt*(s.gammadot+candidate.gammadot);
 s=candidate; count=count+1;
 max_force=max(max_force,abs(force)); minimum_sigma_e=min(minimum_sigma_e,min(cfg.sigma_n-cfg.pa-s.p));
 maximum_coupling_iterations=max(maximum_coupling_iterations,it);
 max_reaction_step=max(max_reaction_step,max(s.reaction_increment));
 minimum_phi=min(minimum_phi,min(s.phi)); maximum_phi=max(maximum_phi,max(s.phi));
 maximum_boundary_p=max(maximum_boundary_p,max(abs(s.p([1,end]))));
 maximum_boundary_T=max(maximum_boundary_T,max(abs(s.T([1,end]))));
 if s.V>peak_V, peak_V=s.V; peak_state=s; end
 if ~event_active && s.V>=0.1, event_active=true; end
 if event_active && s.V<1e-4, event_count=event_count+1; event_active=false; end
 if ~fields_active && ~fields_complete && s.V>=1e-4, fields_active=true; end
 if fields_active && ~fields_complete && s.time-last_field>=field_interval
  if nf>=cfg.field_max
   error('geometric:FieldCapacity','Dense output capacity exceeded; no data thinning permitted.');
  elseif nf==capacity
   capacity=min(cfg.field_max,2*capacity); fn=fieldnames(fields);
   for j=1:numel(fn),fields.(fn{j})(capacity,end)=0;end
  end
  nf=nf+1;
  for j=1:numel(field_names)
   fields.(field_names{j})(nf,:)=single(s.(field_names{j}));
  end
  fields.time_seconds(nf)=s.time; fields.V(nf)=s.V; fields.step_seconds(nf)=dt; fields.tauf(nf)=s.tauf;
  fields.x(nf)=s.x;fields.theta(nf)=s.theta;last_field=s.time;
  if event_count>=1 && isinf(field_stop),field_stop=s.time+cfg.field_tail_seconds;end
  if s.time>=field_stop,fields_complete=true;end
 end
 if mod(count,25)==0 || s.V>=1e-4 || s.time>=cfg.end_seconds
  ih=ih+1;if ih>size(history,1),history(2*size(history,1),18)=0;end
  history(ih,:)=[s.time/31536000,s.V,s.x,s.tauf,max(s.gammadot),max(s.p),dt, ...
    max(s.T),max(s.phi_pd),sum(weights.*s.xi)/h,max(s.xi),min(s.p),min(s.phi),max(s.phi), ...
    min(s.alphahy)/cfg.alpha_hy,max(s.alphahy)/cfg.alpha_hy,force,s.theta];
 end
 elapsed=toc(timer);
 if elapsed-last_progress>=30
  fprintf('PROGRESS case=%s steps=%d days=%.9g V=%.7g events=%d wall_s=%.2f rejected=%d\n', ...
   cfg.case_code,count,s.time/86400,s.V,event_count,elapsed,rejected);
  progress=struct('count',count,'time',s.time,'V',s.V,'event_count',event_count,'wall_seconds',elapsed);
  save(fullfile(out,'progress.mat'),'progress');last_progress=elapsed;
 end
 if fields_complete && cfg.stop_after_first_event,break,end
 dt=min([dt*1.08,cfg.state_step_fraction*cfg.dt_scale*cfg.D_RS/max(s.V,realmin),cfg.dt_limit*cfg.dt_scale]);
end
catch failure
 save(fullfile(out,'failure.mat'),'failure','s','candidate','cfg','dt','count','history','ih');
 rethrow(failure)
end
run_time=toc(timer);completed=s.time>=cfg.end_seconds || (cfg.stop_after_first_event && fields_complete);
history=history(1:ih,:);timev=history(:,1);history_V=history(:,2);history_x=history(:,3);
history_tauf=history(:,4);history_maxgammadot=history(:,5);history_p_max=history(:,6);
history_dt=history(:,7);history_T_max=history(:,8);history_phi_pd_max=history(:,9);
history_xi_mean=history(:,10);history_xi_max=history(:,11);history_p_min=history(:,12);
history_phi_min=history(:,13);history_phi_max=history(:,14);history_alpha_min=history(:,15);
history_alpha_max=history(:,16);history_force_residual=history(:,17);history_theta=history(:,18);
history_f_r=history_tauf./(cfg.sigma_n-cfg.pa-history_p_max);
diagnostics=struct('minimum_sigma_e',minimum_sigma_e,'maximum_force_residual_Pa',max_force, ...
 'peak_V',peak_V,'peak_time',peak_state.time,'event_count',event_count,'fields_complete',fields_complete, ...
 'maximum_boundary_p',maximum_boundary_p,'maximum_boundary_T',maximum_boundary_T, ...
 'accepted_steps',count,'rejected_steps',rejected,'maximum_coupling_iterations',maximum_coupling_iterations, ...
 'maximum_reaction_increment',max_reaction_step,'minimum_phi',minimum_phi,'maximum_phi',maximum_phi);
a=cfg.a;b=cfg.b;V0=cfg.V0;f0=cfg.f0;Lambda=cfg.Lambda;alpha_th=cfg.alpha_th;
alpha_hy=cfg.alpha_hy;rho_c=cfg.rho_c;ny=n;gamma_c=cfg.gamma_c;
gammadotc=cfg.taub0/(cfg.eta*h);fc=f0+(a-b)*log(gammadotc/cfg.gammadot0);
deltac=rho_c*h/(Lambda*fc);tc=cfg.taub0/(cfg.eta*deltac);
accepted_steps=count;final_time=s.time;
save(fullfile(out,'fault_slip_data.mat'),'timev','history_V','history_x','history_tauf', ...
 'history_maxgammadot','history_p_max','history_dt','history_T_max','history_phi_pd_max', ...
 'history_xi_mean','history_xi_max','history_p_min','history_phi_min','history_phi_max', ...
 'history_alpha_min','history_alpha_max','history_force_residual','history_theta','history_f_r', ...
 'diagnostics','parameters','a','b','V0','f0','Lambda','alpha_th','alpha_hy','rho_c','ny','h','yy', ...
 'gamma_c','gammadotc','fc','deltac','tc','run_time','completed','accepted_steps','final_time','-v7');
fn=fieldnames(fields);for j=1:numel(fn),fields.(fn{j})=fields.(fn{j})(1:nf,:);end
spatial=fields;spatial.yy=yy;spatial.h=h;spatial.interval=field_interval;spatial.complete=fields_complete;
spatial.phi_e_initial=initial_state.phi_e;spatial.phi_p_initial=initial_state.phi_p;
save(fullfile(out,'first_event_fields.mat'),'spatial','parameters','-v7.3');
final_state=s;save(fullfile(out,'final_state.mat'),'parameters','completed','final_state','peak_state','initial_state','-v7');
result=struct('completed',completed,'final_days',s.time/86400,'run_time',run_time, ...
 'accepted_steps',count,'samples',ih,'diagnostics',diagnostics);
fprintf('DONE case=%s completed=%d days=%.12g steps=%d wall_s=%.3f peak_V=%.9g events=%d max_force_Pa=%.5g\n', ...
 cfg.case_code,completed,s.time/86400,count,run_time,peak_V,event_count,max_force);
end

function [s,ok,iter,force]=coupled_step(old,dt,c,step_p_T)
s=old;ok=false;force=Inf;
for iter=1:c.max_coupling_iterations
 previous=s;
 if c.decomposition_switch
  integral=.5*dt*c.A*(exp(-c.Q./(c.R*(c.T0+old.T)))+exp(-c.Q./(c.R*(c.T0+s.T))));
  increment=(1-old.xi).*(-expm1(-integral));
 else,increment=zeros(size(old.xi));end
 xi=old.xi+increment;
 if c.dilatancy
  target=c.eps_pl*log(max(s.gammadot,realmin)./c.gammadot_dilatancy_ref);
  phi_pd=old.phi_pd+(target-old.phi_pd).*(-expm1(-dt*s.gammadot/c.gamma_c));
 else,phi_pd=zeros(size(old.phi_pd));end
 plastic=-(phi_pd-old.phi_pd)/(dt*c.Sc);
 [T,p]=step_p_T(old.T,old.p,s.gammadot,s.tauf,dt,s.alphahy,plastic,c.dilatancy,increment/dt);
 if min(c.T0+T)<=0 || min(c.sigma_n-c.pa-p)<=0 || any(~isfinite([T,p])),return,end
 s.T=T;s.p=p;s.xi=xi;s.phi_pd=phi_pd;s.time=old.time+dt;s.reaction_increment=increment;
 s=porosity(s,c);
 if min(s.phi)<=0 || max(s.phi)>=1,return,end
 [s.tauf,s.V,s.gammadot,s.theta,s.x,force]=mechanics(p,old,dt,c,previous.tauf);
 if ~isfinite(s.V) || any(~isfinite(s.gammadot)),return,end
 error_measure=max([max(abs(s.p-previous.p))/1e6,max(abs(s.T-previous.T))/100, ...
    max(abs(s.xi-previous.xi))*100,max(abs(s.phi_pd-previous.phi_pd))*1e4, ...
    max(abs(s.gammadot-previous.gammadot))/(max(s.gammadot)+c.gammadot0), ...
    abs(s.V-previous.V)/(s.V+c.V0)]);
 if error_measure<c.coupling_tolerance && iter>=2
  ok=true;return
 end
end
end

function s=porosity(s,c)
s.phi_e=c.phi0*(c.beta_n*s.p+c.lambda_n*s.T);
s.phi_pc=c.Phi*c.m_d_t*s.xi;
s.phi_p=s.phi_pc+s.phi_pd;s.phi=c.phi0+s.phi_e+s.phi_p;
if c.permeability_evol
 nodal=c.alpha_hy*(s.phi/c.phi0).^3;
 s.alphahy=[nodal(1),2./(1./nodal(1:end-1)+1./nodal(2:end)),nodal(end)];
else,s.alphahy=c.alpha_hy*ones(1,numel(s.p)+1);end
end

function [tau,V,gamma,theta,x,residual]=mechanics(p,old,dt,c,tau_guess)
% Solve the same sinh friction law and quasi-dynamic force balance together.
drive=c.taub0+c.k*(c.loading_rate*(old.time+dt)-old.x-.5*dt*old.V);
eta=c.eta+.5*c.k*dt;
assert(drive>0,'damaged:NegativeDrive','Positive-slip branch lost.');
sigma_e=c.sigma_n-c.pa-p;
% A cellwise rate bound avoids hundreds of Newton iterations when a hot,
% pressurized trial starts from a much stronger previous shear traction.
% theta <= theta_n+dt supplies a conservative upper friction coefficient.
Fmax=c.f0/c.a+(c.b/c.a)*log(c.V0*(old.theta+dt)/c.D_RS);
if strcmpi(c.model,'gaussian')
 L=log(drive/eta)-log(2*c.V0)+Fmax;
 zb=L+log(2);if L<350,zb=asinh(exp(L));end
 bound=c.a*min(sigma_e)*zb;
else
 L=log(drive/eta)-log(2*c.gammadot0*c.weights)+Fmax;
 zb=L+log(2);small=L<350;zb(small)=asinh(exp(L(small)));
 bound=min(c.a*sigma_e.*zb);
end
lo=0;hi=min(drive,bound*(1+1e-12));tau=min(max(tau_guess,lo),hi);
for j=1:80
 V_trial=(drive-tau)/eta;
 midpoint=.5*(old.V+V_trial);zstate=midpoint*dt/c.D_RS;
 if zstate<1e-5
  decay=exp(-zstate);g=1-zstate/2+zstate^2/6-zstate^3/24;
  gp=-.5+zstate/3-zstate^2/8;
 else
  decay=exp(-zstate);g=-expm1(-zstate)/zstate;
  gp=((zstate+1)*decay-1)/zstate^2;
 end
 theta=old.theta*decay+dt*g;
 theta_V=dt/(2*c.D_RS)*(-old.theta*decay+dt*gp);
 F=c.f0/c.a+(c.b/c.a)*log(c.V0*theta/c.D_RS);
 if strcmpi(c.model,'gaussian')
  pressure=max(p);se=c.sigma_n-c.pa-pressure;z=tau/(c.a*se);
  lograte=log(c.V0)+z-F+log(-expm1(-2*z));velocity=exp(lograte);
  dvelocity=velocity*(1/tanh(z)/(c.a*se)+(c.b/c.a)*theta_V/(eta*theta));
  gamma=velocity*c.gaussian_shape;
 else
  z=tau./(c.a*sigma_e);
  lograte=log(c.gammadot0)+z-F+log(-expm1(-2*z));gamma=exp(lograte);
  velocity=sum(c.weights.*gamma);
  derivative=gamma.*(1./tanh(z)./(c.a*sigma_e)+(c.b/c.a)*theta_V/(eta*theta));
  dvelocity=sum(c.weights.*derivative);
 end
 residual=tau+eta*velocity-drive;
 if abs(residual)<=c.force_tolerance*max(drive,1)
  V=velocity;x=old.x+.5*dt*(old.V+V);
  % theta was evaluated at force-balance V_trial; their discrepancy is
  % bounded by force_tolerance*drive/eta and is monitored in validation.
  return
 end
 if residual>0,hi=tau;else,lo=tau;end
 next=tau-residual/(1+eta*dvelocity);
 if ~isfinite(next) || next<=lo || next>=hi,next=.5*(lo+hi);end
 tau=next;
end
error('damaged:MechanicalSolve','Force balance did not converge.');
end

function c=defaults(c)
input_names=fieldnames(c);
for i=1:numel(input_names)
 if isnumeric(c.(input_names{i})),c.(input_names{i})=double(c.(input_names{i}));end
end
aliases={'phi_0','phi0';'beta','Sc';'Lambda_n','lambda_n';'Dc','D_RS';'rect_m_frac','mbar';'ending_time','end_seconds';'Ta','T0';'V_L','loading_rate'};
for i=1:size(aliases,1)
 if isfield(c,aliases{i,1}) && ~isfield(c,aliases{i,2}),c.(aliases{i,2})=c.(aliases{i,1});end
end

required={'output_directory','case_code','phi0','beta_n','lambda_n','Sc','Lambda','alpha_hy', ...
 'alpha_th','rho_c','sigma_n','pa','T0','Phi','m_d_t','Pr','Er','mbar'};
for i=1:numel(required),assert(isfield(c,required{i}),'Missing cfg.%s',required{i});end
d=struct('model','mrsf','width',26e-6,'ny',300,'h',.0005,'end_seconds',60*86400, ...
 'dt_scale',1,'dt_limit',51.84,'max_steps',Inf,'stop_after_first_event',false,'random_seed',334500, ...
 'a',.01,'b',.015,'f0',.6,'V0',1e-6,'D_RS',.005,'gamma_c',1,'eps_pl',1.7e-4, ...
 'eta',4.5e6,'R_k',.5,'loading_rate',1e-6,'A',10^15.47,'Q',319e3,'R',8.3145, ...
 'boundary_type','dirichlet','state_step_fraction',1/32,'pressure_step',7e4,'temperature_step',5, ...
 'reaction_step',.01,'log_velocity_step',.08,'max_coupling_iterations',30, ...
 'coupling_tolerance',2e-6,'force_tolerance',2e-12,'field_interval',.001,'field_max',50000, ...
 'field_tail_seconds',.05,'first_event_dt_max',1e-5,'dense_velocity_threshold',.01);
names=fieldnames(d);for i=1:numel(names),if ~isfield(c,names{i}),c.(names{i})=d.(names{i});end,end
c.case_code=char(c.case_code);assert(ismember(c.case_code,{'0','000','001','101','011','111'}),'Invalid case code');
if strcmpi(c.model,'gaussian')
 y=linspace(-c.h/2,c.h/2,c.ny);shape=exp(-y.^2/(2*c.width^2))/(sqrt(2*pi)*c.width);
 % The resolved shape integrates to V on the finite gouge grid.
 c.gaussian_shape=shape/trapz(y,shape);
end
end

function report=selftest(cfg)
% Numerical identities independent of a production trajectory.
cfg.output_directory=tempdir;cfg.case_code='111';c=defaults(cfg);
n=41;c.ny=n;c.yy=linspace(-c.h/2,c.h/2,n);dy=c.h/(n-1);
c.weights=dy*ones(1,n);c.weights([1,end])=dy/2;c.gammadot0=c.V0/c.h;
c.taub0=c.f0*(c.sigma_n-c.pa);c.k=c.R_k*(c.sigma_n-c.pa)*(c.b-c.a)/c.D_RS;
c.permeability_evol=true;
s=struct('p',linspace(0,40e6,n),'T',linspace(0,700,n), ...
 'xi',linspace(0,1,n),'phi_pd',zeros(1,n));s=porosity(s,c);
assert(all(s.phi_e>=0));assert(max(abs(s.phi_pc-c.Phi*c.m_d_t*s.xi))<1e-14);
assert(all(s.alphahy>=c.alpha_hy));
nodal=c.alpha_hy*(s.phi/c.phi0).^3;
assert(max(abs(s.alphahy(2:n)-2./(1./nodal(1:n-1)+1./nodal(2:n))))<1e-18);
old=struct('time',1e5,'x',.1,'theta',c.D_RS/c.V0,'V',c.V0);
force_error=0;state_error=0;
for pressure=[0 30e6 60e6]
 for dt=[0 1e-4 1]
  [tau,V,rate,theta,x,res]=mechanics(pressure*ones(1,n),old,dt,c,c.taub0);
  equilibrium=c.taub0+c.k*(c.loading_rate*(old.time+dt)-x)-tau-c.eta*V;
  force_error=max(force_error,abs(equilibrium));
  assert(abs(V-sum(c.weights.*rate))<1e-10*max(V,c.V0));
  assert(abs(res)<1e-3 && abs(equilibrium)<1e-3);
  vavg=.5*(old.V+V);z=vavg*dt/c.D_RS;
  expected=old.theta*exp(-z)+dt*(-expm1(-z))/max(z,realmin);
  if dt==0,expected=old.theta;end
  state_error=max(state_error,abs(theta-expected)/max(expected,1));
  assert(abs(theta-expected)<1e-6*max(expected,1));
 end
end
% With sealed boundaries, uniform sources obey the integrated balances.
robin=robin_boundary_config(c.alpha_hy,c.alpha_th,c.rho_c,dy);
robin.type='robin';
robin.hydraulic_diffusivity=0;robin.thermal_conductivity=0;
fixed=struct('ny',n,'alpha_th',c.alpha_th,'rho_c',c.rho_c,'Lambda',c.Lambda, ...
 'dy',dy,'rect_m_frac',c.mbar,'Pr',c.Pr,'Er',c.Er,'robin',robin, ...
 'temperature_ambient',0,'pressure_ambient',0);
step=make_p_T_stepper(fixed);dt=.01;increment=.002;dphi=.00001;
Tn=400*ones(1,n);pn=2e6*ones(1,n);rate=.2*ones(1,n);tau=20e6;
[T,p]=step(Tn,pn,rate,tau,dt,c.alpha_hy,-dphi/(c.Sc*dt),1,increment/dt);
expectedT=Tn+dt*tau*rate/c.rho_c-c.mbar*c.Er*increment;
expectedp=pn+c.Lambda*(T-Tn)-dphi/c.Sc+c.mbar*c.Pr*increment;
thermal_error=max(abs(T-expectedT));pressure_error=max(abs(p-expectedp));
assert(thermal_error<1e-6 && pressure_error<.01);
report=struct('passed',true,'force_error_Pa',force_error,'relative_state_error',state_error, ...
 'sealed_thermal_error_K',thermal_error,'sealed_pressure_error_Pa',pressure_error);
disp(report)
end
