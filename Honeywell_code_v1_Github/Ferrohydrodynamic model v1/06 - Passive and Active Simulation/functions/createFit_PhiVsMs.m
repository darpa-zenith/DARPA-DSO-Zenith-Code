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
[xData, yData]  = prepareCurveData( Ms_data, phi_data ); %Prepare data 
% inputs for curve fitting

% Set up fittype and options.

% defining the fit type, linear
ft              = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
% additional curvefitting options
opts            = fitoptions( 'Method', 'NonlinearLeastSquares' );          
opts.Display    = 'Off';
opts.StartPoint = 0.0002141;

% Fit model to data.
[fitresult, gof]= fit( xData, yData, ft, opts );




