function Pi_m = gridinterp_Pi_m(x_m,z_m,H_m,params)
% Computes the gridded interpolant for the magnetic potential. 
%
% Pi_m = gridinterp_Pi_m(x_m,z_m,H_m,params)
%
% Inputs:
%     1st Input  : x_m
%     2nd Input  : z_m
%     3rd Input  : H_m
%     4th Input  : params
%
% Outputs:
%     Pi_m  
%
% V1.1 - Álvaro Romero-Calvo (01/01/2024)

% Set vector of H values
H       = linspace(0,1e6,10000);
Pi_m_v  = zeros(length(H),1);

% Integrate
for i = 2:length(H)
    Pi_m_v(i) = -params.mu0/params.rho * trapz(H(1:i),magnetization(H(1:i),params));
end

% Gridded interpolant
Pi_m_H  = griddedInterpolant(H,Pi_m_v,'linear');
Pi_m    = griddedInterpolant(x_m,z_m,Pi_m_H(H_m),'linear');