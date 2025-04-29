function PI = calc_PI(Z,H_norm)
%CALC_PI Calculate the potential of the body force (magnetic force, and gravity) of the ferrofluid assuming the magnetization of the ferrofluid follows a Langevin curve 
%   Input:    Z - Height [m]
%                H_norm - Strength of magnetic field [A/m] 
%       
%   Output: Total potential of the body force

% Ferrofluid property
chi0         = 0.6733;                        % Initial susceptibility
Ms           = 5607;                     % Saturation magnetization    [A/m]
rho           = 1060.9;                    % Density  [kg/m^3]      
gamma    = 3*chi0/Ms;            % Langevin curve parameter

% Calculate the potential of the body forces
PI_m       =  -(4*pi*1e-7)/rho*Ms* (abs(gamma * H_norm) - log(abs(gamma * H_norm)) + log(1 - exp(-2 * abs(gamma * H_norm))) - log(2))/ gamma;
PI_g        =  9.806*Z;
PI            =  PI_g+PI_m;
end