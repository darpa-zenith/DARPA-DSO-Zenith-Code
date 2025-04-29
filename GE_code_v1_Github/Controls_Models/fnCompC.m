function C = fnCompC(H0, Hmnc, Hmns)

% The equations in this function are based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

% Function to compute system's C matrix, Pg 136 after Eqn(59), [Azhar 2008]

tempHmnc = Hmnc';
tempHmnc = tempHmnc(:);
tempHmns = Hmns';
tempHmns = tempHmns(:);

tempHmncs = [tempHmnc, tempHmns];
tempHmncs = tempHmncs';
tempHmncs = tempHmncs(:);

tempH0mncs = [H0'; tempHmncs];

tempC = [tempH0mncs 0*tempH0mncs];
tempC = tempC';

C = tempC(:)';