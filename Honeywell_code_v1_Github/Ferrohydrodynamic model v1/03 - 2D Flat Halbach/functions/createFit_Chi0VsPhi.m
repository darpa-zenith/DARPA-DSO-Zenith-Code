function [fitresult, gof] = createFit_Chi0VsPhi(phi_data, Chi0_data)
% This code generates a curvefit for the case of Chi0 vs Phi, using a list of
% chi0 and phi data
%
% [fitresult, gof] = createFit_Chi0VsPhi(phi_data, Chi0_data)
%
% Data for fit:
%     1st Input  : phi_data
%     2nd Input  : Chi0_data
%
% Output:
%     fitresult  : a fit object representing the fit.
%     gof        : structure with goodness-of fit info.
%
% V1.1 - Tianyang Hu(08/01/2024)


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( phi_data, Chi0_data );

% Set up fittype and options.
ft             = fittype( 'a*x^b', 'independent', 'x', 'dependent', 'y' );
excludedPoints = xData > 5;
opts           = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display   = 'Off';
opts.StartPoint= [0.5147 1.469];
opts.Exclude   = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData, excludedPoints );
% legend( h, 'Chi0_data vs. phi_data', 'Excluded Chi0_data vs. phi_data', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'phi_data', 'Interpreter', 'none' );
% ylabel( 'Chi0_data', 'Interpreter', 'none' );
% grid on


