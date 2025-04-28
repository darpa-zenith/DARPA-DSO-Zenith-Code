function Y = fnCompY(lambda, z, params)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute Y, refer Eqn (35), [Azhar 2008]

mu = params.mu;
mu0 = params.mu0;
chi = params.chi;

[alpha, beta] = fnCompAlphaBeta(lambda, params);

Y = -(mu/mu0 * beta./alpha + chi./alpha) .* cosh(lambda*z) + ...
    (mu/mu0 * (alpha./beta - chi./alpha) - chi^2./(alpha.*beta)) .* sinh(lambda*z);