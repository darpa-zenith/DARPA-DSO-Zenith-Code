function [Pi_m, dPi_m_dx, dPi_m_dz] = magpot(x,z,params)
% Computes the magnetic potential and its derivatives at selected points
% r,z
%
% [Pi_m, dPi_m_dx, dPi_m_dz] = magpot(x,z,params)
%
% Inputs:
%     1st Input  : x
%     2nd Input  : z
%     3rd Input  : params
%
% Outputs:
%     1st Output : Pi_m 
%     2nd Output : dPi_m_dx
%     3rd Output : dPi_m_dz
%
% V1.1 - Álvaro Romero-Calvo (01/01/2024)

% Parameters
eps_   = 1e-10;
Pi_m_gi= params.Pi_m;

% Compute results
x        = x';
z        = z';
Pi_m     = Pi_m_gi(x,z);
dPi_m_dx = (Pi_m_gi(x+eps_,z)-Pi_m_gi(x-eps_,z))/2/eps_;
dPi_m_dz = (Pi_m_gi(x,z+eps_)-Pi_m_gi(x,z-eps_))/2/eps_;

% Convert to meshgrid
Pi_m     = Pi_m';
dPi_m_dx = dPi_m_dx';
dPi_m_dz = dPi_m_dz';