function [T,p]=legacy_robin_p_T(Tn,pn,gamma,tauf,dt,Tnorth,Tsouth, ...
    pnorth,psouth,M_alpha_hy,alpha_th,rho_c,Lambda,dy,plastic,use_plastic, ...
    reaction,mbar,Pr,Er,boundary)
%LEGACY_ROBIN_P_T Extend archived interior equations with Robin edge rows.
% Keep the archived source/face indexing, so changing the boundary option
% does not silently replace the production interior discretization.
% The current THMC* solver uses the collocated make_p_T_stepper instead.
n=numel(Tn);
[ht,hp]=robin_transfer_coefficients(boundary,rho_c);
gamma=nodal(gamma,n);reaction=nodal(reaction,n);plastic=nodal(plastic,n);
mask=nodal(double(use_plastic),n);
heat=tauf*gamma/rho_c-mbar*Er*reaction;
heat(2:n-1)=tauf*gamma(1:n-2)/rho_c-mbar*Er*reaction(1:n-2);
T=solve_robin_diffusion(Tn,alpha_th,heat,dt,dy,ht,[Tsouth Tnorth]);
source=Lambda*(T(:)-Tn(:))/dt+mask.*plastic+mbar*Pr*reaction;
dT=T(:)-Tn(:);
source(2:n-1)=Lambda*dT(2:n-1)/dt ...
    +mask(2:n-1).*plastic(1:n-2)+mbar*Pr*reaction(2:n-1);
if isscalar(M_alpha_hy),faces=M_alpha_hy;else,faces=M_alpha_hy(1:n-1);end
p=solve_robin_diffusion(pn,faces,source,dt,dy,hp,[psouth pnorth]);
end

function x=nodal(x,n)
if isscalar(x),x=repmat(x,n,1);end
assert(numel(x)==n,'Expected a scalar or one value per physical node.');
x=x(:);
end
