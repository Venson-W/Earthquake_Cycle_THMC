function gammadot = shear_strain(tauf, theta, gammadot0, V0, f0, b, a, Dc, sigma_n, p)
% GAMMADOT Computes the strain rate from shear stress using RSF inverse law.
%
% Inputs:
%   tauf        - shear stress (scalar or array)
%   theta       - state variable (same size as tauf)
%   gammadot0   - reference strain rate (scalar)
%   f0          - reference friction coefficient (scalar)
%   b, a        - rate-and-state friction parameters (scalars)
%   Dc          - critical slip distance (scalar)
%   sigma_n     - normal stress (scalar or array)
%   p           - pore pressure (scalar or array)
%
% Output:
%   strainrate  - strain rate (same size as tauf)

    effective_stress = sigma_n - p;

    argument = 2 .* gammadot0 .* sinh(tauf ./ (a .* effective_stress)) ...
               .* exp(-f0 ./ a - b./a * log(V0*theta./Dc));

    gammadot = argument;
end
