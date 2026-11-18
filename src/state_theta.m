function theta = state_theta(thetan, gammadot, dt, Dc_prime)
%STATE_THETA Computes the state variable theta based on aging law exact solution
%
% Inputs:
%   thetan - initial state variable theta (scalar or matrix)
%   gammadot     - shear strain rate (matrix)
%   dt     - time step (scalar)
%   Dc_primer    - critical slip distance (scalar)
%
% Output:
%   theta  - updated state variable after time dt

    theta = thetan .* exp(-gammadot .* dt ./ Dc_prime) + Dc_prime ./ gammadot .* (1 - exp(-gammadot .* dt ./ Dc_prime));
end
