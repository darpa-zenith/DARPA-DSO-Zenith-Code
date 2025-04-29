function [fitresult, gof] = createFit_PhiVsMs(Ms_data, phi_data)
% This code generates a curvefit for the case of Phi vs Ms, using a list of
% Ms and phi data
%
% [fitresult, gof] = createFit_PhiVsMs(Ms_data, phi_data)
%
% Data for fit:
%     1st Input  : Ms_data
%     2nd Input  : phi_data
%
% Output:
%     fitresult  : a fit object representing the fit.
%     gof        : structure with goodness-of fit info.
%
% V1.1 - Tianyang Hu(08/01/2024)


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( Ms_data, phi_data );

% Set up fittype and options.
ft             = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
opts           = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display   = 'Off';
opts.StartPoint= 0.0002141;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'phi_data vs. Ms_data', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'Ms_data', 'Interpreter', 'none' );
% ylabel( 'phi_data', 'Interpreter', 'none' );
% grid on


