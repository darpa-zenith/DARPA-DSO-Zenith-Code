%% This is the main script to create state space model, transfer function 
% and step response for MFDM (Magnetic Ferrofluid Deformable Mirror) system

% The model is developed based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

clear all; close all; clc;


%%
% Truncation parameters for Bessel function of first kind
M_bessel = 4; % an integer value >= 0
N_bessel = 4; % an integer value > 0

% Location of deformation measurement in polar coordinations
r_in        = 0;
theta_in    = 0;

% Location of EM coil in polar coordinates
r_coils_in        = 0;
theta_coils_in    = 0;

% Roots of Bessel function of first kind, J_m(E) = 0, m= 0,1,2, ...

E_orig = [2.405,  5.520,  8.654, 11.792, 14.931, 18.071, 21.212, 24.352, 27.493, 30.635;
          3.832,  7.016, 10.173, 13.324, 16.471, 19.616, 22.760, 25.904, 29.047, 32.19;
          5.136,  8.417, 11.620, 14.796, 17.960, 21.117, 24.270, 30.569, 33.717, 36.863; 
          6.380,  9.761, 13.015, 16.223, 19.409, 22.583, 25.748, 28.908, 32.065, 35.219;
          7.588, 11.065, 14.373, 17.616, 20.827, 24.019, 27.199, 30.371, 33.537, 36.699; 
          8.772, 12.339, 12.339, 15.700, 18.980, 22.218, 25.430, 28.627, 31.812, 34.989;
          9.936, 13.589, 17.004, 20.321, 23.586, 26.820, 30.034, 33.233, 36.422, 39.603;
         11.086, 14.821, 18.288, 21.642, 24.935, 28.191, 31.423, 34.637, 41.031, 44.215;
         12.225, 16.038, 19.555, 22.945, 26.267, 29.546, 32.796, 36.026, 39.240, 42.444;
         13.354, 17.241, 20.807, 24.234, 27.584, 30.885, 34.154, 37.400, 40.629, 43.844];

E = E_orig(1:M_bessel+1, 1:N_bessel);
% Load parameters and compute lambda and z
run('modelParameters.m');

R = MdlParams.R;
lambda = E./R;  % Eqn 43, [Azhar 2008]

h = MdlParams.EMcoil2filmDist;
z = -h;
% Compute Y and Z functions

Y = fnCompY(lambda, z, MdlParams); % function to compute Y, refer Eqn (35), [Azhar 2008]
Z = fnCompZ(lambda, z, MdlParams); % function to compute Z, refer Eqn (36), [Azhar 2008]
% Compute omega^2, F, Hmnc, Hmns, omega_d

omega2                              = fnCompOmega2(Y, Z, lambda, MdlParams); % function to compute omega^2, refer Eqn (50), [Azhar 2008]
Fmn                                 = fnCompF(lambda, Y, E, MdlParams); % function to compute Fmn, refer Eqn (52), [Azhar 2008]
[H0_coils, Hmnc_coils, Hmns_coils]  = fnCompH_coils(lambda, E, r_coils_in, theta_coils_in); % function to compute H_mnc^j and H_mns^j,
                                                                                            % refer Eqn (53)  and Eqn(56), [Azhar 2008]
omega_d                             = fnCompOmegad(lambda, MdlParams); % function to compute omega_d_mn, refer Eqn (58), [Azhar 2008]
[H0, Hmnc, Hmns]                    = fnCompH(lambda, r_in, theta_in); % function to compute H_mnc and H_mns, refer Eqn (44), [Azhar 2008]
% Compute A, B, C, D matrices

A = fnCompA(omega2, omega_d); % function to compute system's A matrix, Pg 135 after Eqn(59), [Azhar 2008]
B = fnCompB(Fmn, H0_coils, Hmnc_coils, Hmns_coils, length(r_coils_in)); % function to compute system's B matrix, Pg 135 & Pg 136after Eqn(59), [Azhar 2008]
C = fnCompC(H0, Hmnc, Hmns); % function to compute system's C matrix, Pg 136 after Eqn(59), [Azhar 2008]
D = zeros(length(r_in), length(r_coils_in)); % D matrix is 0 matrix
% Compute Open Loop Transfer function, Bode plot and Step response

G = tf(ss(A,B,C,D));

%% Ad-hoc Controller to improve the system response - Tuning
kp = 100;
ki = 0.3;

K = tf([kp, ki],[1 0]);
H = K*G/(1+K*G);

% figure(3)
% bode(H)
options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.figure(1)

subplot(121);
bode(G,{0.1 100},options)
% xlim([0 300])
grid on
title('OL Bode Diagram')
subplot(122);
bode(H,{0.1 100},options)
% xlim([0 300])
title('CL Bode Diagram')
grid on;


%% Simulink model

out = sim("MFDM_Mdl_SISO.slx");

% Simulink post process
out.logsout; % read log data from simulink

Open_loop = get(out.logsout,"OL");
Open_loop_response = Open_loop.Values;

% without saturator
Close_loop1 = get(out.logsout,"CL_1");
Close_loop_response1 = Close_loop1.Values;


figure(4)
plot(Open_loop_response,'Linewidth',1.5)
hold on;
plot(Close_loop_response1,'Linewidth',1.5)
grid on
title('Step Response of FDM system');
xlabel('Time (s)');
ylabel('Normalized Amplitude');
legend('OL',['kp_s=',num2str(kp),' ki_s=',num2str(ki)]);

