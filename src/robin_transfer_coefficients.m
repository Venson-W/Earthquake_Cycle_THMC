function [temperature_transfer, pressure_transfer] = robin_transfer_coefficients(robin, rho_c)
%ROBIN_TRANSFER_COEFFICIENTS Normalized Robin exchange factors [m/s].
%   Temperature flux is thermal_conductivity/layer_thickness times the
%   boundary temperature change from ambient [W/m^2]. Dividing by rho_c
%   gives the flux entering the temperature evolution equation.
%   The pressure equation uses hydraulic_diffusivity/layer_thickness.
%   Equal host/gouge reference storage is assumed, as in the paper.
%   Output ordering is [south north].
%   type='dirichlet' returns infinite transfer for both fields so the
%   diffusion solver imposes the supplied ambient values exactly.
%   An absent type field uses the default Dirichlet behavior.

validateattributes(robin, {'struct'}, {'scalar'});
validateattributes(rho_c, {'numeric'}, {'real','finite','scalar','positive'});
boundary_type = 'dirichlet';
if isfield(robin, 'type')
    boundary_type = robin.type;
end
if ~(ischar(boundary_type) && isrow(boundary_type)) && ...
        ~(isstring(boundary_type) && isscalar(boundary_type))
    error('robin:InvalidType', 'robin.type must be ''robin'' or ''dirichlet''.');
end
if ~any(strcmpi(boundary_type, {'robin','dirichlet'}))
    error('robin:InvalidType', 'robin.type must be ''robin'' or ''dirichlet''.');
end
if strcmpi(boundary_type, 'dirichlet')
    temperature_transfer = [Inf Inf];
    pressure_transfer = [Inf Inf];
    return
end
names = {'thermal_conductivity','hydraulic_diffusivity', ...
    'thermal_layer_thickness','hydraulic_layer_thickness'};
values = cell(1,4);
for k = 1:4
    if ~isfield(robin, names{k})
        error('robin:MissingSetting', 'Missing robin.%s.', names{k});
    end
    value = robin.(names{k});
    validateattributes(value, {'numeric'}, {'real','finite','vector','nonnegative'});
    if ~ismember(numel(value), [1 2]) || (k >= 3 && any(value <= 0))
        error('robin:InvalidSetting', ...
            'robin.%s must be a scalar or [south north] pair; thicknesses must be positive.', names{k});
    end
    if isscalar(value), value = [value value]; end
    values{k} = reshape(value,1,2);
end
temperature_transfer = values{1} ./ (rho_c .* values{3});
pressure_transfer = values{2} ./ values{4};
end
