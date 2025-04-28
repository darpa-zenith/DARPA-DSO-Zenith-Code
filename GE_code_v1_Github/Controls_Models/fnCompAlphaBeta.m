function [alpha, beta] = fnCompAlphaBeta(lambda, params)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute alpha and beta, refer Eqn (37) and Eqn (38), [Azhar 2008]

mu  = params.mu;
mu0 = params.mu0;
d   = params.thickness;

alpha = tanh(lambda*d) - coth(lambda*d);

beta = (mu/mu0)*tanh(lambda*d) - coth(lambda*d);