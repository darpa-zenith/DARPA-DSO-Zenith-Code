function [fitresult, gof] = createFit_RhoVsPhi(phi_data, rho_data)
% This code generates a curvefit for the case of rho vs Phi, using a list of
% rho and phi data
%
% [fitresult, gof] = createFit_RhoVsPhi(phi_data, rho_data)
%
% Data for fit:
%     1st Input  : phi_data
%     2nd Input  : rho_data
%
% Output:
%     fitresult  : a fit object representing the fit.
%     gof        : structure with goodness-of fit info.
%
% V1.1 - Tianyang Hu(08/01/2024)
%% Fit: 'untitled fit 1'.
[xData, yData]      = prepareCurveData( phi_data, rho_data );%Prepare data
% inputs for curve fitting

% Set up fittype and options.

% defining the fit type, linear
ft                  = fittype( 'a*x+1000', 'independent', 'x', 'dependent', 'y' ); 
% additional curvefitting options
opts                = fitoptions( 'Method', 'NonlinearLeastSquares' );              
opts.Display        = 'Off';
opts.StartPoint     = 0.392227019534168;

% Fit model to data.
[fitresult, gof]    = fit( xData, yData, ft, opts );



