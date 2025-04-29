function omega_d        = fnCompOmegad(lambda, params)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute omega_d_mn, refer Eqn (58), [Azhar 2008]

rho     = params.rho;
eta     = params.eta;

omega_d    = 4 * eta/rho * lambda.^2;