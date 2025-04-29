function res = residuals_setPI(Y,params)
% Residual computation function for the determination of the 2D 
% equipotential interface shape. 
%
% res = residuals_interface(Y,params)
% 
% Inputs:
%           Y       Y coordinate
%           params  Parameters structure
%
% Outputs:
%           res     residual, defined as Pi_iterated - Pi_target
% V1.1 - Álvaro Romero-Calvo (01/01/2024)

% Expand parameters structure
x_m         = params.x_ff;  
% y_m   = params.y_m;
Pi          = params.Pi;
eqp_height  = params.delta/2;

% Compute potential at desired point
Pi_0  = Pi(0,eqp_height);

% Evaluate
res   = Pi(x_m,Y) - Pi_0;