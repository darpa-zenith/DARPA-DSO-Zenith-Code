function F = fnCompF(lambda, Y, E, params)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute Fmn, refer Eqn (52), [Azhar 2008]

d       = params.thickness;
h       = params.EMcoil2filmDist;
R       = params.R;
rho     = params.rho;
chi     = params.chi;
B0      = params.HHcoilB0;

[M, N] = size(E);

tempBesselTerm = zeros(size(E));
for ii=1:M    
    if ii==1
        kVal = 1;
    else
        kVal = 2;
    end
    tempBesselTerm(ii,:) = kVal./besselj(ii, E(ii,:)).^2;
end

F = -(chi/rho) * B0 * tanh(lambda*d)./Y .*lambda.^2 .* (1/(pi*R^2)) .* tempBesselTerm;