function [H0, Hc, Hs]    = fnCompH_coils(lambda, E, r, theta)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute H_mnc^j and H_mns^j, refer Eqn (53)  and Eqn(56), [Azhar 2008]


H0 = besselj(0, kron(r, lambda(1,:))) .* cos(kron(theta, zeros(size(lambda(1,:)))));


[M, N] = size(E);

Hc = zeros(M-1, N*length(r));
Hs = zeros(M-1, N*length(r));
for ii=2:M
    Hc(ii-1, :) = besselj(ii-1, kron(r, lambda(ii,:))) .* cos(kron(theta, (ii-1)*ones(size(lambda(1,:)))));
    Hs(ii-1, :) = besselj(ii-1, kron(r, lambda(ii,:))) .* sin(kron(theta, (ii-1)*ones(size(lambda(1,:)))));
end

