%% Error vs lambda and h/lambda plots:
%  This file plots the force, xi, and Rosenzwieg instability vs magnet wavelength and ferrofluid height-to-wavelength ratio.
%  v2.4, Eric Comstock, 1230 ET 05/04/2023

clear
clc
clearvars
close all% these four lines close all existing plots and free up memory and variables to be re-allocated later

%% Code configuration
N                  = 101;                            % defining number of points per axis - This code should perfom at around 12 seconds per myriad points
addpath('functions');                                % defining functionpath
%% Parameters
params.sigma       = 0.025;                          % Surface tension (N/m), for water based ferrofluid
params.mu0         = 1.25663706212e-6;               % Vacuum permeability
params.M0          = 1.48/params.mu0;                % N42 (1.32) ->N52 (1.48) grade magnet, gauss is in kg/(A*s^2), this is A/m
params.delta       = 0.001;                          % Thickness of the ferrofluid layer (m)
params.x           = 0;                              % here x (m) and y (m) represent the simulation point
params.y           = params.delta;                   % y (m) is defined as the distance above the bottom of the ferrofluid layer
params.boverlam    = 1;                              % Magnet thickness b over magnet wavelength lambda - held constant here (unitless)
params.g_target    = 10;                             % target number of g's that we want our solutions to have - this is done my varying the Ms of the ferrofluid

%% Constants
params.g           = 9.81;                           % standard gravity of Earth (m/s^2)
params.mu0         = 4*pi*1e-7;                      % vacuum permeability (m*kg/(s^2 A^2))

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

%% Follow-up parameters
lam_1D    = 0.0254 * 4 * 2 .^ linspace(-8.1, 4.1, N);% 1D list of magnet wavelengths
hr_1D     = linspace(0, 1, N);                       % 1D ratio of ferrofluid height to magnet wavelength
[lam, hr] = meshgrid(lam_1D, hr_1D);                 % 2D input variables

% Because I do a Newton method to find Ms, I have to use for loops. error, Ms, xi, and keps all start out as zero matrices to be popuated later.
error     = zeros(N);                                % Surface error
Ms        = zeros(N);                                % Saturation magnetization - needed for Rosenswieg constraint
xi        = zeros(N);                                % Calculating xi for verification - must be much less than 1 for our approximations to work
keps      = zeros(N);                                % Calculating k*epsilon for verification - must be much less than 1 for our approximations to work
inst      = zeros(N);                                % Calculating difference between magnetization and critical magnetization

for i = 1:N
    for j = 1:N
        [error(i, j), Ms(i, j), xi(i, j), keps(i, j), inst(i, j)] = find_E(lam(i, j), hr(i, j),RhovsMs_cfit, Chi0vsRho_cfit, params);   % finding error and optimum Ms
        progress = 'progress at ' + string(100*(i * N - N + j) / N^2) + '%';                                % progress report on screen
        progress                                                                                            % printing the progress in % for ease of determining how much time is left
    end
end

manuf     = lam - 0.0254 /16 * 4;                    % This is positive when the magnets can be manufactured now, and negative when new techniques are needed

%% Generating plots

figure(1)                                                                       % generating figure to apply everything to
hold on                                                                         % applying hold so that multiple plots can be stacked
megacontour(log10(lam_1D / 0.0254), hr_1D, inst, 0, -1000, '#1b9e77');       % Drawing Rosenswieg constraint
megacontour(log10(lam_1D / 0.0254), hr_1D, xi, 0.1, 0.085, '#7f7f7f');          % drawing xi constraint
megacontour(log10(lam_1D / 0.0254), hr_1D, hr, 0.1, 0.115, 'k');                % drawing h/lambda constraint
leg=legend('Rosenswieg instability constraint', '', 'ξ constraint', '', 'h/λ constraint', '');% initializing legend
clines = reshape([1; 3] * 10.^ linspace(-10, 9, 20), [1, 40]);                  % contour line formation - these are the places where lines will be drawn to allow presentation of results in decimal notation

[c11, c12] = contour(log10(lam / 0.0254), hr, error * 1e6, clines,'ShowText','on','Color', '#7570b3');% WFE contours
c12.DisplayName = 'Mid-wavefront error (μm)';                                   % adding this contour to the legend
c12.LineWidth = 1.5;                                                            % making this contour readable
clabel(c11, c12, 'Color', [0.459 0.439 0.702]);                                 % coloring text

[c21, c22] = contour(log10(lam / 0.0254), hr, 1./(lam / 4 / 0.0254), 2 .^ linspace(-10,10,21),'ShowText','on','Color', 'k');%
c22.DisplayName = 'Magnets per inch';                                           % adding this contour to the legend
c22.LineWidth = 1.5;                                                            % making this contour readable

[c31, c32] = contour(log10(lam / 0.0254), hr, Ms, clines,'ShowText','on','Color', '#1b9e77');
c32.DisplayName = 'Saturation Magnetization (A/m)';                             % adding this contour to the legend
c32.LineWidth = 1.5;                                                            % making this contour readable
clabel(c31, c32, 'Color', [0.106 0.620 0.467]);                                 % coloring text

[c41, c42] = contour(log10(lam / 0.0254), hr, xi, 10 .^ [-3 -2 -1 0 2 4 6 8],'ShowText','on','Color', '#d95f02');
c42.DisplayName = 'ξ';                                                          % adding this contour to the legend
c42.LineWidth = 1.5;                                                            % making this contour readable
clabel(c41, c42, 'Color', [0.851 0.373 0.008]);                                 % coloring text
% clabel(c41, c42, 'Color', [0.906 0.161 0.541]);%'#e7298a' %this is a different color from colorblender than can be applied to new contour plots

% This code block is plotting kε - which is extremely low (<10^-9) and thus not useful to plot most of the time. Uncomment if needed.
% [c51, c52] = contour(log10(lam), hr, log10(keps),'ShowText','on','Color', '#66a61e');
% c52.DisplayName = 'log_{10}(kε)';                 % adding this contour to the legend
% c52.LineWidth = 1.5;                              % making this contour readable
% clabel(c51, c52, 'Color', [0.400 0.651 0.118]);   % coloring text

[c61, c62] =contourf(log10(lam / 0.0254), hr, 0-manuf, [0 0], "FaceAlpha",0.15);% paint all points that cannot be manufactured
leg.String(end) = [];                                                           % remove this from legend
colormap("autumn");                                                             % paint them red
t=text(-2,hr_1D(N)*1.05,{'Requires new manufacturing'},'r');                    % stating that these magnets need new manufacturing
t.Color = [1 0 0];                                                              % color the text red
t.FontSize = 14;                                                                % make the text bigger

xlabel("log_{10}(Magnet wavelength λ (in))")                                    % labelling x-axis
ylabel('Ferrofluid height h/λ (unitless)')                                      % labelling y-axis
title('log_{10}(Mid-wavefront error)')                                          % adding title
hold off                                                                        % un-applying hold to retrun to normal - that way I can make more plots in the future if I want to

%% Computing specific values - uncomment if nessesary

% '1/8 in magnet, h/λ = 1/8'
% [gm, em] = find_E(0.0254 / 8 * 4, 1/8,RhovsMs_cfit);
% gm * 1e9
% em
% 
% '1/16 in magnet, h/λ = 1/8'
% [gm, em] = find_E(0.0254 / 16 * 4, 1/8,RhovsMs_cfit);
% gm * 1e9
% em
% 
% '1/32 in magnet, h/λ = 1/8'
% [gm, em] = find_E(0.0254 / 32 * 4, 1/8,RhovsMs_cfit);
% gm * 1e9
% em

%% Plot format
set(findobj('Type', 'Legend'), 'Interpreter', 'tex', 'box','on');
set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
    'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'tex','TickDir','out');
set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'tex');
set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [2,2,25,15]);
set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'tex');