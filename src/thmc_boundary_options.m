function boundary=thmc_boundary_options(options,alpha_hy,alpha_th,rho_c,h)
%THMC_BOUNDARY_OPTIONS Dirichlet by default; Robin is explicitly selected.
boundary=robin_boundary_config(alpha_hy,alpha_th,rho_c,h/299);
boundary.type='dirichlet';
if isfield(options,'boundary_type')
    boundary.type=validatestring(options.boundary_type,{'dirichlet','robin'});
end
if isfield(options,'boundary_layer_thickness')
    boundary.hydraulic_layer_thickness=options.boundary_layer_thickness;
    boundary.thermal_layer_thickness=options.boundary_layer_thickness;
end
if isfield(options,'robin')
    validateattributes(options.robin,{'struct'},{'scalar'});
    names=fieldnames(options.robin);
    allowed={'hydraulic_diffusivity','thermal_conductivity', ...
        'hydraulic_layer_thickness','thermal_layer_thickness'};
    for i=1:numel(names)
        assert(ismember(names{i},allowed),'Unknown Robin option: %s',names{i});
        boundary.(names{i})=options.robin.(names{i});
    end
end
robin_transfer_coefficients(boundary,rho_c);
end
