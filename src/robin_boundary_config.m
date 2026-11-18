function robin = robin_boundary_config(alpha_hy, alpha_th, rho_c, reference_spacing)
%ROBIN_BOUNDARY_CONFIG Fixed host properties for heat and fluid exchange.
%   alpha_hy: initial gouge hydraulic diffusivity [m^2/s].
%   alpha_th: initial gouge thermal diffusivity [m^2/s].
%   rho_c: gouge volumetric heat capacity [J/(m^3 K)].
%   reference_spacing: physical boundary-layer thickness [m].
%
% type is 'dirichlet' by default; 'dirichlet' fixes both boundary values to
% the supplied ambient changes. Each physical field may subsequently be
% set to a scalar or [south north] pair.
% Thermal conductivity is rho_c*alpha_th, not alpha_th itself.
% Keep these properties fixed when the gouge permeability evolves.
% The drivers use h/(300-1), the original grid interval, as the default
% physical layer thickness. Hold it fixed during grid refinement.

validateattributes(alpha_hy, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(alpha_th, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(rho_c, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(reference_spacing, {'numeric'}, {'real','finite','scalar','positive'});
robin.type = 'dirichlet';
robin.hydraulic_diffusivity = alpha_hy;
robin.thermal_conductivity = rho_c * alpha_th;
robin.hydraulic_layer_thickness = reference_spacing;
robin.thermal_layer_thickness = reference_spacing;
end
