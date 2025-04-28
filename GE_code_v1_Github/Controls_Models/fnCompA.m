function A = fnCompA(omega2, omega_d)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute system's A matrix, Pg 135 after Eqn(59), [Azhar 2008]

[M, N] = size(omega2);

A0 = zeros(2*N);
for ii=1:N
    A0(2*(ii-1)+1:2*(ii-1)+2, 2*(ii-1)+1:2*(ii-1)+2) = [0 1; -omega2(1,ii) -omega_d(1,ii)];
end

Acs = zeros(4*N*(M-1));
for ii=2:M
    for jj=1:N
        atomA = [0 1; -omega2(ii,jj) -omega_d(ii,jj)];
        atomAcs = blkdiag(atomA, atomA);
        Acs(4*N*(ii-2)+(4*(jj-1)+1:4*(jj-1)+4), 4*N*(ii-2)+(4*(jj-1)+1:4*(jj-1)+4)) = atomAcs;
    end
end

A = blkdiag(A0, Acs);
