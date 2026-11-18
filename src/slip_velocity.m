function V = slip_velocity(gammadot, yy_gammadot)
% SLIP_VELOCITY Computes slip velocity by integrating strain rate over space.
%
% Inputs:
%   gammadot     - strain rate values (array)
%   yy_gammadot  - spatial coordinates corresponding to gammadot (array)
%
% Output:
%   V            - slip velocity (scalar), computed via spatial integration

    %V = trapz(yy_gammadot, gammadot, 2);
    V = trapz(yy_gammadot, gammadot, 2);
end
