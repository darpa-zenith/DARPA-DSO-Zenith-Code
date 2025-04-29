function [H0, Hmnc, Hmns] = fnCompH(lambda, r, theta)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute H_mnc and H_mns, refer Eqn (44), [Azhar 2008]

H0 = besselj(0, lambda(1,:)*r) * cos(0*theta);

[M, N] = size(lambda);
Hmnc = zeros(M-1, N);
Hmns = zeros(M-1, N);
for ii=2:M
    Hmnc(ii-1, :) = besselj(ii-1, lambda(ii,:)*r) * cos((ii-1)*theta);
    Hmns(ii-1, :) = besselj(ii-1, lambda(ii,:)*r) * sin((ii-1)*theta);
end

% for ii=2:M
%     Hmnc(ii-1, :) = besselj(ii-1, lambda(ii,:)*r) * cos(0*theta);
%     Hmns(ii-1, :) = besselj(ii-1, lambda(ii,:)*r) * sin(0*theta);
% end