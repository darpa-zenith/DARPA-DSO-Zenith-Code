function [fitresult, gof] = createFit_PhiVsMs(Ms_data, phi_data)
%CREATEFIT(MS_DATA,PHI_DATA)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: Ms_data
%      Y Output: phi_data
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 24-Dec-2023 12:17:57


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( Ms_data, phi_data );

% Set up fittype and options.
ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.0002141;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );



