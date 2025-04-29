%% Roughness vs error plots
%  This file plots the force, xi, and Rosenzwieg instability vs ferrofluid saturation magnetization and ferrofluid height.
%  It plots multiple magnets, so the plots may be compared.
%  v5.1, Eric Comstock, 1230 ET 05/04/2024 AD

clear
clc
clear
close all%these four lines close all existing plots and free up memory and variables to be re-allocated later

%% Code Configuration
stack_plots = 1; % If true, stacks all of the plots into one. If false, breaks them into four.
N           = 1000;                           % Number of points per axis
addpath('functions');                         % Adding function path
%% Constants
mu0         = 4*pi*1e-7;                      % Permeability of free space (m*kg/(s^2 A^2))
g0          = 9.81;                           % Gravity acceleration (m/s^2)

%% Magnetic Curve Interpolation
% Interpolation data for EMG water-based series ferrofluid (Ferrotec Inc.)
% URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/
% EMG   304, 308 408 507  508  509 605  607  700   705  707  708

% Data
Chi0_data   = [5.03 0.5  0.5  1.63 0.88 0.5  3.02 1.63 12.57 4.04 1.51 0.63];                % Initial Susceptibility
Ms_data     = [27.5 6.6  6.6  11   6.6  3.3  22   11   35.5  22   11   6.6] * 795.7747;      % Saturation Magnetization [A/m]
phi_data    = [4.5  1.2  1.2  2    1.2  0.6  3.9  2    5.8   3.9  2    1.2];                 % Magnetic Particle concentration [%]
rho_data    = [1240 1060 1070 1120 1070 1030 1180 1100 1290  1190 1100 1080];                % Density @ 25 C [kg/m^3]

% Find fit curves
[RhovsMs_cfit, gof1]   = createFit_RhovsMs(Ms_data, rho_data);
[Chi0vsRho_cfit, gof2] = createFit_Chi0vsRho(rho_data, Chi0_data);

%% Parameters
sigma       = 0.025;                                                            % Surface tension (N/m), for water based ferrofluid
M0          = 1.48/mu0;                                                         % Magnet magnetization (A/m), N42 (1.32) ->N52 (1.48) grade magnet, gauss is in kg/(A*s^2), this is A/m
delta       = 0.001;                                                            % Thickness of the ferrofluid layer (m)
width_list  = 1/2*0.0254; [1/2, 1/2, 1, 4, 8, 1, 4, 8, 1, 4, 8] * 0.0254 / 8;   % Magnet widths (m)
b_list      = 1/2*0.0254; [1/2, 8, 1, 1, 1, 4, 4, 4, 8, 8, 8] * 0.0254 / 8;     % Magnet heights (m)
Ms_max      = 37000;                                                            % Maximum saturation magnetization (A/m)
Ms_min      = 100;                                                              % Minimum saturation magnetization (A/m)
h_max       = 0.07;                                                             % Maximum ferrofluid separation (m)
h_min       = 0;                                                                % Minimum ferrofluid separation (m)
sep_mag     = 0.025*0.0254; 2.5e-4;                                             % Magnet separation (m)

%% Follow-up parameters
Ms_1D       = 10.^linspace(log10(Ms_min), log10(Ms_max), N);   % Range of Ms values we want
h_1D        = linspace(h_min,h_max,N);                         % Range of h values - ferrofluid height (m) - we want
[Ms, h]     = meshgrid(Ms_1D, h_1D);                           % 2D grid of Ms and h for use in further calculations
rho         = reshape(RhovsMs_cfit(Ms),size(Ms));              % kg/m^3; interpolated earlier and in function createFit_RhovsMs
chi0        = reshape(Chi0vsRho_cfit(rho),size(Ms));           % interpolated earlier and in function createFit_Chi0vsRho

%% Compute variables for magnet of choice. 
t = cell(length(b_list));
for i = 1:length(b_list)                      % Generating figure titles
    t{i}=string(width_list(i) / 0.0254) + " x " + string(b_list(i) / 0.0254) + " in magnet";
end

for i = 1:length(b_list)                      % Different figures for each magnet
    width       = width_list(i);              % Width of magnet
    b           = b_list(i);                  % Height of magnet - both this and the width are from the width and height lists

    % Intermediate results
    lambda      = 4*(width + sep_mag);                                      % Frequency length (m)
    k           = 2*pi./lambda;                                             % Wave number (1/m)
    H0          = (2^1.5 / pi) * M0 * exp(-k*h) * (1 - exp(-k*b));          % Calculating magnetic H field at the surface of the ferrofluid
    b51      = (1 - exp(-5 * k * b))/(10 - 10 * exp(-1 * k * b));           % beta51 - added from Dr. Gomez' calculations
    b91      = (1 - exp(-9 * k * b))/(9 - 9 * exp(-1 * k * b));             % beta91 - added from the new calculations
    xi          = Ms.* exp(k*delta)./H0;                                    % Calculating xi for verification
    
    % Surface force
    f_vol_y     = -mu0*k*Ms.*H0.*exp(-k.*delta) .* (1 + xi);                % Magnetic force in the y direction (N/m3)
    f_vol_y_g   = f_vol_y./(rho*g0);                                        % Magnetic force per unit mass - in units of g
        
    % Surface mid-WFE roughness
    Omega       = (4*sigma*k + rho.*g0./k)./(mu0*Ms.*H0.*exp(-k*delta));    % Quotient in surface error term 
    gam_s_g_M_H = 1 + Omega - 10 * b51 .* exp(-4 * k * (h + delta)) + (9 * b91 - 2 * xi .* b51 .^2) .* exp(-8 * k * (h + delta));%capital gamma in the new expressions
    epsilon     = abs(xi./(4*k*gam_s_g_M_H) .* (1 - 12 * b51 .* (exp(k * delta) - 1) .* exp(-4 * k * (h + delta)) + (b51 .^ 2 + 14 * (exp(k * delta) - 1) * (b91 - 3 * b51 .^ 2)) .* exp(-8 * k * (h + delta))   )); % Mid-wavefront error (m)
    
    % Rosensweig instability
    Hx      = H0;                                                                                               % Worst case scenario for Mallinson configuration - H0 is entirely in the normal direction
    gamma   = 3*chi0./Ms;                                                                                      % Gefinition of gamma
    r0      = sqrt(((H0+Ms.*(coth(gamma.*H0)-1./(gamma.*H0))).*(1+Ms.*(1./(H0.^2.*gamma)-gamma.*csch(H0.*gamma).^2)))./H0);%(geometric mean of the chord and tangent permeabilities)
    Mc      = sqrt(2./mu0.*(1+1./r0).*(sqrt(g0*rho.*sigma)+mu0*(1-r0).^2/2./(1+r0).*Hx.^2));                % Critical magnetization[A/m]
    M       = Ms.*(coth(gamma.*H0)-1./gamma./H0);                                                               % Magnetization of ferrofluid [A/m]
    inst    = M - Mc;                                                                                           % Represents the instability - if this is negative, we are good

    allcon  = max((-10 - f_vol_y_g),0) .* max((1 - xi * 10),0) .* max((0 - inst), 0) .* max((h - 0.1*width), 0);% This is all of our constraints except the Roughness, and is >0 when we are feasible
    
    % Plotting
    figure(i)                                                                                                   % generating the figure for this particular magnet configuration
    hold on
    if stack_plots
        title(t(i))
        megacontour(log10(Ms_1D), h_1D * 1000, f_vol_y_g, -10, -10.5,'#1b9e77');                                % force constraint
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.01, 0.0095,'#823901');                                     % stringent xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.1, 0.095,'#d95f02');                                       % weak xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, inst, 0, -250,'#909090');                                        % Rosensweig constraint
        megacontour(log10(Ms_1D), h_1D * 1000, h, 0.1*width, 0.1*width + 0.02*max(h_1D),'#000000');             % Height constraint

        contourf(log10(Ms_1D), h_1D * 1000, allcon ./ allcon + allcon, [0 0], 'FaceAlpha',0.25)                 % paint all feasible points a translucent color
        colormap("summer")                                                                                      % paint them green
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        legend('Force > 10g constraint', '', '\xi < 0.01 constraint', '', '\xi < 0.1 constraint', '', 'Rosensweig instability constraint', '', 'Height constraint')
        
        clines = reshape([1; 5] * 10.^ linspace(-10, 4, 15), [1, 30]);                                          % contour lines

        [c11, c12] = contour(log10(Ms), h * 1000, epsilon * 1e6, clines,'Color', '#7570b3','ShowText','on');    % plot of log Roughness
        c12.DisplayName = 'Mid-wavefront error (μm)';                                                           % adding this contour to the legend
        c12.LineWidth = 1.5;                                                                                    % making the lines visible
        clabel(c11, c12, 'Color', [0.459 0.439 0.702]);                                                         % coloring the labels

        [c21, c22] = contour(log10(Ms), h * 1000, 0-f_vol_y_g, clines,'Color', '#1b9e77','ShowText','on');      % plot of log force
        c22.DisplayName = 'Magnetic force (g)';                                                                 % adding this contour to the legend
        c22.LineWidth = 1.5;                                                                                    % making the lines visible
        clabel(c21, c22, 'Color', [0.106 0.620 0.467]);                                                         % coloring the labels 

        [c31, c32] = contour(log10(Ms), h * 1000, xi, clines,'Color', '#d95f02','ShowText','on');               % plot of log xi
        c32.DisplayName = '\xi';                                                                                % adding this contour to the legend
        c32.LineWidth = 1.5;                                                                                    % making the lines visible
        clabel(c31, c32, 'Color', [0.851 0.373 0.008]);                                                         % coloring the labels
    else
        sgtitle(t(i))                                                                                           % showing the viewer which magnet this figure is for
    
        ax4 = subplot(2,2,4);                                                                                   % plot of all of the constraints - ax4 is for imposing a different colormap
        bpoint = megacontour(log10(Ms_1D), h_1D * 1000, f_vol_y_g, -10, -10.5,'#1b9e77');                       % force constraint and ideal point
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.01, 0.0095,'#823901');                                     % stringent xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.1, 0.095,'#d95f02');                                       % weak xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, inst, 0, -250,'#404040');                                        % Rosensweig constraint
    
        contourf(log10(Ms_1D), h_1D * 1000, allcon ./ allcon + allcon, [0 0], "FaceAlpha",0.25)                 % paint all feasible points a translucent color
        colormap(ax4, "summer")                                                                                 % paint them green
    
        axis([min(log10(Ms_1D)) max(log10(Ms_1D)) min(h_1D * 1000) max(h_1D * 1000)])                           % making sure the axes in this subplot matches those of the others
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        legend('Force > 10g constraint', '', '\xi < 0.01 constraint', '', '\xi < 0.1 constraint', '', 'Rosensweig instability constraint')
        title('Plot of constraints')
    
        subplot(2,2,1)          % plot of log Roughness
        contourf(log10(Ms), h * 1000, log10(epsilon),'ShowText','on');
        colorbar;
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        title('log_{10}(Mid-wavefront error (m))')
    
        subplot(2,2,2)          % plot of log force
        contourf(log10(Ms), h * 1000, log10(0-f_vol_y_g),'ShowText','on');
        colorbar;
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        title('log_{10}(Magnetic force (g))')
    
        subplot(2,2,3)          % plot of log xi
        contourf(log10(Ms), h * 1000, log10(xi),'ShowText','on');
        colorbar;
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        title('log_{10}(\xi)')
    end
    hold off
    drawnow
end

%% Plot format
set(findobj('Type', 'Legend'), 'Interpreter', 'tex', 'box','on');
set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
    'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'tex','TickDir','out');
set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'tex');
set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [2,2,25,15]);
set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'tex');