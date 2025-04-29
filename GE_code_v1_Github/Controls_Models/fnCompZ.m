function Z = fnCompZ(lambda, z, params)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute Z, refer Eqn (36), [Azhar 2008]

mu = params.mu;
chi = params.chi;

[alpha, beta] = fnCompAlphaBeta(lambda, params);

Z = (beta.*cosh(lambda*z) + chi*sinh(lambda*z)) * chi/mu .* 1./alpha;
% Z = (beta.*cosh(lambda*z) + chi*sinh(lambda*z)) * chi/1.78 .* 1./alpha;
