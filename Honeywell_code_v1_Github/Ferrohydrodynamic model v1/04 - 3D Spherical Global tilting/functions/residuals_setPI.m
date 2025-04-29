function res = residuals_setPI(Z, z0, params)
% Residual computation function for the determination of the 3D 
% equipotential interface shape. Each residual vector element is one (x,y)
% point of the interface
%
% res = residuals_interface(Z, z0, params)
% 
% Inputs:
%           Z       Vector of state (elevations of (x,y) points)
%           z0      Target height at (0,0)
%           params  Parameters structure
%
% Outputs:
%           res     Vector of residuals Pi_iterated - Pi_target

% Expand parameters structure
x_m   = params.x_m;  
y_m   = params.y_m;
Pi    = params.Pi;

% Evaluation grid for ferrofluid surface
x_t       = x_m(:,:,1);
y_t       = y_m(:,:,1);

% Compute potential at desired point
Pi_0  = Pi(0,0,z0);

% Evaluate
res   = Pi(x_t,y_t,Z) - Pi_0;