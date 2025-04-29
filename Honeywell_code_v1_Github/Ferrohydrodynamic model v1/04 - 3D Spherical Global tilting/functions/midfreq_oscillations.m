function epsilon = midfreq_oscillations(x,delta,params)
% Computes the surface roughness produced by the Mallinson Configuration +
% first two harmonics of Halbach array assuming h/lambda > 0.1 and xi << 1
%
% epsilon = roughness(x,f,params)
%
% Inputs:
%     1st Input  : x
%     2nd Input  : delta
%     3rd Input  : params
%
% Outputs:
%     epsilon    : epsilon
%
% V1.1 - Álvaro Romero-Calvo (06/27/2024)

% Parameters
M0     = params.M0;
h      = params.h;
k      = params.k;
b      = params.b;
Ms     = params.Ms;
rho    = params.rho;
g      = params.g;
sigma  = params.sigma;
mu0    = params.mu0;
DNoN   = params.DNoN;
chi0   = params.chi0;

% H0
H0     = 2*sqrt(2)/pi * M0 * exp(-k*h) * (1 - exp(-k*b));

% xi
xi     = Ms/H0 * exp(k*delta);

% Omega
Omega  = (4.*sigma.*k + rho.*g./k)./(mu0.*Ms.*H0.*exp(-k.*delta));

% Gamma
Gamma  = 1./(1 + Omega + DNoN/2 * (1 + xi/2 * exp(-k*delta) * (1+DNoN/2)));

% F1 and F2
F1     = 2/chi0; 
F2     = 2/chi0;

% epsilon
epsilon= Gamma./k .* (DNoN/2 + xi/4 .* (1 - DNoN*F1...
    - 1/2 * DNoN^2 * (1 + F2 - exp(-k*delta))));

