function [PI] = Pi(B, params)
    %% Pi calculates the magnetic potential from magnetic flux density B
    %
    % Inputs: 
    % B:            Magnetic Flux Density Norm                  [T]
    % params        Parameter variable structure
    %
    % Outputs: 
    % PI            Magnetic mass-specific potential function   [J/kg]
    %
    % V4.1, Eric Comstock, 01/04/2024 AD
    
    %% Parameters
    mu0     = params.mu0;                                   % Vaccum Pemeability [N/A^2]
    rho     = params.rho;                                   % Ferrofluid material density [kg/m^3]                     
    Ms      = params.Ms;                                    % Saturation Magnetization of ferrofluid [A/m]
    gamma   = params.gamma;                                 % Ferrofluid-specific property [m/A]
    
    %% Calculating Magnetic Field H
    H_norm  = 1./mu0.*B;                                      % Magnetic field [A/m] (M = [0;0] assumed)	
    input   = H_norm .* gamma;                               % Input to ln(sinh) calculation
    ls_out  = abs(input) - log(abs(input)) + log(1 - exp(-2 .* abs(input))) - log(2); % ln(sinh) calculation
    
    %% Calculate Mass Force Potential PI
    PI      = -mu0./rho.*Ms.*(ls_out/gamma);                   %Magnetic potential through the Langevin formula
end