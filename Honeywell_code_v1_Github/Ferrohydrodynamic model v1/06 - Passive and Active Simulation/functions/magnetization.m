function M = magnetization(H,params)
% Returns the magnetization vector module in A/m for a given external 
% magnetic field H in A/m
%
% M = magnetization(H,params)
%
% Inputs:
%     1st Input  : H
%     2nd Input  : params
%
% Outputs:
%     M  
%
% V1.1 - Álvaro Romero-Calvo (01/01/2024)

% Parameters
chi0     = params.chi0;
Ms       = params.Ms;

% Langevin Curve
gamma    = 3*chi0/Ms;
M        = Ms*(coth(gamma*H)-1/gamma./H);
M(H==0)  = 0;
end