function theta = state_theta(thetan, V, dt, Dc)
%STATE_THETA Computes the state variable theta based on aging law exact solution
%
% Inputs:
%   thetan - initial state variable theta (scalar or matrix)
%   V     - shear strain rate (matrix)
%   dt     - time step (scalar)
%   Dcr    - critical slip distance (scalar)
%
% Output:
%   theta  - updated state variable after time dt

    theta = thetan .* exp(-V .* dt ./ Dc) + Dc ./ V .* (1 - exp(-V .* dt ./ Dc));
end
