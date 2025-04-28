function B = fnCompB(Fmn, H0, Hcmn, Hsmn, Ncoils)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute system's B matrix, Pg 135 and Pg 136 after Eqn(59), [Azhar 2008]

[M, N] = size(Fmn);

F0 = zeros(2*N, N);
for ii=1:N
    F0(2*(ii-1)+(1:2), ii) = [0; Fmn(1, ii)];
end

Fcs = zeros(4*N*(M-1), 2*N*(M-1));
for ii=2:M
    for jj=1:N
        atomF = [0; Fmn(ii,jj)];
        atomFcs = blkdiag(atomF, atomF);
%         Fcs(4*(ii-2)+1:4*(ii-2)+4, 2*(ii-2)+1:2*(ii-2)+2) = atomFcs; 
        Fcs(4*N*(ii-2)+(4*(jj-1)+1:4*(jj-1)+4), 2*N*(ii-2)+(2*(jj-1)+1:2*(jj-1)+2)) = atomFcs; 
    end
end
F = blkdiag(F0, Fcs);

H0      = reshape(H0', numel(H0)/Ncoils, Ncoils);
Hcsmn   = zeros(2*N*(M-1), Ncoils);
for ii=1:Ncoils
    tempHcmn = Hcmn(:, (ii-1)*N+(1:N));
    tempHsmn = Hsmn(:, (ii-1)*N+(1:N));

    tempHcmn = tempHcmn';
    tempHcmn = tempHcmn(:);
    tempHsmn = tempHsmn';
    tempHsmn = tempHsmn(:);

    tempHcsmn = [tempHcmn, tempHsmn]';
    tempHcsmn = tempHcsmn(:);

    Hcsmn(:,ii)  = tempHcsmn;
end
H = [H0; Hcsmn];


B = F * H;