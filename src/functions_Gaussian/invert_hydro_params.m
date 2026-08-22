function [beta_n, Lambda_n, beta_f, Lambda_f, eta_f, k] = invert_hydro_params(n, Lambda, beta, alpha_hy)
%INVERT_HYDRO_PARAMS  Invert bulk hydro-thermo parameters for porous rock.
%
%   Inputs:
%       n        - porosity (scalar or array)
%       Lambda   - bulk thermal pressurization coefficient (Λ)
%       beta     - bulk compressibility (β) of the fluid–rock system [Pa^-1]
%       alpha_hy - hydraulic diffusivity α_hy [m^2/s] (or consistent units)
%
%   Constants (fluid properties):
%       beta_f   = 4.5e-10 Pa^-1
%       Lambda_f = 3e-4   1/K
%       eta_f    = 1e-3   Pa·s   (water viscosity)
%
%   Relations used:
%       β        = n * (β_f + β_n)
%       Λ        = (Λ_f - Λ_n) / (β_f + β_n)
%       α_hy     = k / (η_f * n * (β_f + β_n))
%     which implies:
%       β_n      = β / n - β_f
%       Λ_n      = Λ_f - Λ * (β_f + β_n)
%       k        = α_hy * β * η_f
%
%   Outputs:
%       beta_n   - rock (skeleton) compressibility [Pa^-1]
%       Lambda_n - rock (skeleton) thermal expansivity [1/K]
%       k        - permeability [m^2] (if α_hy in m^2/s, β in Pa^-1, η_f in Pa·s)

    % --- Constants ---
    beta_f   = 4.5e-10;   % Pa^-1
    Lambda_f = 3e-4;      % 1/K
    eta_f    = 1e-3;      % Pa·s

    % --- Rock-frame compressibility β_n ---
    beta_n = beta ./ n - beta_f;

    % --- Rock-frame thermal expansivity Λ_n ---
    Lambda_n = Lambda_f - Lambda .* (beta_f + beta_n);

    % --- Permeability k ---
    % Using α_hy = k / (η_f * n * (β_f + β_n)) and β = n*(β_f+β_n) ⇒ k = α_hy * β * η_f
    k = alpha_hy .* beta .* eta_f;
end
