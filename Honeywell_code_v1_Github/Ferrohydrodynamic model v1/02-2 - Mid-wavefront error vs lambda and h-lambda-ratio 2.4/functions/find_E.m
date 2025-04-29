function [e, Ms, xi, keps, inst] = find_E(lambda, hratio,RhovsMs_cfit, Chi0vsRho_cfit, params)
%FIND_E This function the WFE in meters, the Minimum saturation magnetizartion of the ferrofluid needed to reach 10g, the xi parameter, and k*epsilon for a particular magnet configuration
% Thus, by finding the Minimum saturation magnetizartion of the ferrofluid needed to reach 10g, it can calculate the optimum for the other inputs of the
% find_GandE function. It is used by plot_E to actually generate the plots.

% bpoint = megacontour(x, y, Z, c0, c1, color)
%
% Inputs:   Any parameters we might want to change
% lambda            Wavelength of magnets
% hratio            height of ferrofluid divided by lambda
% RhovsMs_cfit      Fitting array for the relation between density and Ms for diluted EMG-700
% params            list of parameters needed for the code to function - inputted in structure format
%
% Outputs:
% Ms                Minimum saturation magnetizartion of the ferrofluid needed to reach 10g
% e                 Mid-wavefront error in meters
% xi                Xi parameter        - must be <<1 for our approximations to work
% keps              k*epsilon parameter - must be <<1 for our approximations to work
% inst              Becomes positive when the Rosensweig instability is active
%
% V2.2, Eric Comstock, 1500 ET 12/2/2023 AD
    Ms                  = 28500;                    % initial guess: EMG700 ferrofluid
    g_target            = params.g_target;          % Target force to achieve by varying the ferrofluid magnetization
    
    dMs                 = 0.001;                    % this is somewhat arbitrary, so I just made it relatively small

    %% Newton method
    [g, e, xi, keps]    = find_GandE(lambda, hratio, Ms,RhovsMs_cfit, Chi0vsRho_cfit, params);
    i = 0;                                          % counts iterations, and will cut the loop if it goes over 500
    while (abs(g + g_target) > 1e-9 && i < 500)     % main Newton method loop
        [g_high, e, xi, keps, inst] = find_GandE(lambda, hratio, Ms + dMs,RhovsMs_cfit, Chi0vsRho_cfit, params);  % finding parameters when Ms is slightly higher
        [g_low, e, xi, keps, inst] = find_GandE(lambda, hratio, Ms - dMs,RhovsMs_cfit, Chi0vsRho_cfit, params);   % finding parameters when Ms is slightly lower
        dg = (g_high - g_low) / (2 * dMs);          % derivative of force in g with respect to Ms
        Ms = Ms - (g + g_target) / dg;              % Next Newton method step
        [g, e, xi, keps, inst]   = find_GandE(lambda, hratio, Ms,RhovsMs_cfit, Chi0vsRho_cfit, params);           % re-finding center point for derivative evaluation after a step
        i = i+1;                                    % intreasing iteration number by 1
    end

    if i>500                                        % this whole part is to avoid an error and return data that obviously shows that an error is occuring
        g = 0;
        e = 0;
        xi = 0;
        keps = 0;
        inst = 1;
    end
end