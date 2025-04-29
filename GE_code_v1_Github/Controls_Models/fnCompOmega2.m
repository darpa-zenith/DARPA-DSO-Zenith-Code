function omega2 = fnCompOmega2(Y,Z,lambda,params)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute omega^2, refer Eqn (50), [Azhar 2008]

d       = params.thickness;
g       = params.g;
sigma   = params.sigma;
rho     = params.rho;
chi     = params.chi;
B0      = params.HHcoilB0;


omega2  = g * tanh(lambda*d) .* lambda +...
          sigma/rho * tanh(lambda*d) .* lambda.^3 + ...
          chi/rho * B0^2 * tanh(lambda*d) .* lambda.^2 .* Z./Y;