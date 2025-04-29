function [f_vol_y_g, Roughness, xi, keps, inst] = find_GandE(lambda, hratio, Ms,RhovsMs_cfit, Chi0vsRho_cfit, params)
%FIND_GANDE This function finds the magnetic force in g's, the WFE in meters, the xi parameter, and k*epsilon for a particular magnet configuration

% bpoint = megacontour(x, y, Z, c0, c1, color)
%
% Inputs:   Any parameters we might want to change
% lambda            Wavelength of magnets
% hratio            height of ferrofluid divided by lambda
% Ms                Saturation magnetization of ferrofluid
% RhovsMs_cfit      Fitting array for the relation between density and Ms for diluted EMG-700
% params            list of parameters needed for the code to function - inputted in array format and unpacked with the unpack function
%
% Outputs:
% f_vol_y_g         Magnetic force on the ferrofluid in g's
% Roughness         Mid-wavefront error in meters
% xi                Xi parameter        - must be <<1 for our approximations to work
% keps              k*epsilon parameter - must be <<1 for our approximations to work
% inst              Becomes positive when the Rosensweig instability is active
%
% V2.2.1, Eric Comstock, 1700 ET 02/04/2024 AD
    
    %% Constants and Parameters
    g           = params.g;                                             % Gravity acceleration (m/s^2)
    mu0         = params.mu0;                                           % Permeability of free space (m*kg/(s^2 A^2))
    sigma       = params.sigma;                                         % Surface tension (N/m), for water based ferrofluid
    M0          = params.M0;                                            % Magnet magnetization (A/m), N42 (1.32) ->N52 (1.48) grade magnet, gauss is in kg/(A*s^2), this is A/m
    delta       = params.delta;                                         % Thickness of the ferrofluid layer (m)
    x           = params.x;                                             % here x (m) and y (m) represent the simulation point
    y           = params.y;                                             % y (m) is defined as the distance above the bottom of the ferrofluid layer
    boverlam    = params.boverlam;
    
    rho         = RhovsMs_cfit(Ms);                                     % interpolated earlier and in function createFit_RhovsMs
    chi0        = Chi0vsRho_cfit(rho);                                  % interpolated earlier and in function createFit_Chi0vsRho

    %% list of our four magnets' widths and heights (b) - modify this to add/remove magnets
    
    b           = lambda * boverlam;                                    % Height of magnet - from finding the correct magnet thickness from lambda and a nondimensional parameter from params
    h           = hratio*lambda;                                        % range of h values - ferrofluid height (m) - we want
    
    %% Surface roughness and force
    
    k           = 2.*pi./lambda;                                        % Wave number (1/m)
    H0          = (2^1.5 / pi) * M0 .* exp(-k*h) .* (1 - exp(-k*b));    % Calculating magnetic H field at the surface of the ferrofluid
    b51      = (1 - exp(-5 * k * b))/(10 - 10 * exp(-1 * k * b));       % beta51 - added from Dr. Gomez' calculations
    b91      = (1 - exp(-9 * k * b))/(9 - 9 * exp(-1 * k * b));         % beta91 - added from the new calculations
    xi          = Ms.* exp(k.*delta)./H0;                               % Calculating xi for verification
    f_vol_term1 = -mu0.*k.*Ms.*H0.*exp(-k.*y);                          % first term for magnetic force, which is broken into two terms to simplify the calculation
    f_vol_term2 = (1 + xi .* exp(2.*k.*(y-delta)).*cos(2*k.*x));        % second term for ''       ''  - see above
    f_vol_y     = f_vol_term1 .* f_vol_term2;                           % Magnetic force in the y direction (N/m3)
    
    f_vol_y_g   = f_vol_y./(rho*g);                                     % Magnetic force per unit mass - in units of g
    
    Omega       = (4.*sigma.*k + rho*g./k)./(mu0*Ms.*H0*exp(-k*delta)); % Quotient in surface error term 
    gam_s_g_M_H = 1 + Omega - 10 * b51 .* exp(-4 * k * (h + delta)) + (9 .* b91 - 2 .* xi .* b51 .^2) .* exp(-8 .* k .* (h + delta));%capital gamma in the new expressions
    Roughness   = abs(xi./(4*k*gam_s_g_M_H) .* (1 - 12 * b51 * (exp(k * delta) - 1) .* exp(-4 * k * (h + delta)) + (b51 .^ 2 + 14 * (exp(k * delta) - 1) * (b91 - 3 * b51 .^ 2)) .* exp(-8 * k * (h + delta))   )); % Mid-wavefront error (m)
    keps        = k* Roughness;                                         % k*eps for verification (unitless)

    %% Rosenswieg instability
    Hx      = H0;                                                                                               % Worst case scenario for Mallinson configuration - H0 is entirely in the normal direction
    gamma   = 3*chi0./Ms;                                                                                      % Gefinition of gamma
    r0      = sqrt(((H0+Ms.*(coth(gamma.*H0)-1./(gamma.*H0))).*(1+Ms.*(1./(H0.^2.*gamma)-gamma.*csch(H0.*gamma).^2)))./H0);%(geometric mean of the chord and tangent permeabilities)
    Mc      = sqrt(2./mu0.*(1+1./r0).*(sqrt(g*rho.*sigma)+mu0*(1-r0).^2/2./(1+r0).*Hx.^2));                % Critical magnetization[A/m]
    M       = Ms.*(coth(gamma.*H0)-1./gamma./H0);                                                               % Magnetization of ferrofluid [A/m]
    inst    = M - Mc;                                                                                           % Represents the instability - if this is negative, we are good

end