%% This is the main script to create state space model, transfer function 
% and step response for MFDM (Magnetic Ferrofluid Deformable Mirror) system

% The model is developed based on the work of Dr. Azhar Iqbal's 2008 paper 
% titled "Modeling and Exerimental Evaluation of a Circular Magnetic-Fluid 
% Deformable Mirror" in International Journal of Optomechatronics. 
% This paper is referred to as [Azhar 2008] in comments below

clear all; 
close all; 
clc;

%%
% Truncation parameters for Bessel function of first kind
M_bessel = 4; % an integer value >= 0
N_bessel = 4; % an integer value > 0

% Load parameters 
run('modelParameters.m');

%% User Defined selections

userDefinedCoilNums   = [1, 2, 3, 4, 5, 6, 7];

%% Location of EM coil in polar coordinates

coil_loc_radius = 35;
No_of_coils = 7;
r_coils_in_all        = [0,repmat(coil_loc_radius,1,No_of_coils-1)] * mm2m;
% Theta locations 
min_angle = 0;
max_angle = 300;
theta_coils_in_all    = [0,linspace(min_angle,max_angle,No_of_coils-1)]*deg2rad;

r_coils_in_all        = r_coils_in_all(userDefinedCoilNums);
theta_coils_in_all    = theta_coils_in_all(userDefinedCoilNums);


% assuming all displacement location are above coils
r_in_all = r_coils_in_all;
theta_in_all = theta_coils_in_all;
%% Raidus of influence of active coils
Min_distance = coil_loc_radius * mm2m;

%% Plot surface displacement measurement and coils layout
coil_Location = [r_coils_in_all',theta_coils_in_all'];
% figure(1)
for i = 1:1:length(coil_Location)
    [x,y] = pol2cart(coil_Location(i,2),coil_Location(i,1));
    % r = 0.01;
    % circle1(x,y,r);
    % text(x,y,num2str(userDefinedCoilNums(i)));
    % grid on
end
% hold on
Location_disp = [r_in_all',theta_in_all'];
for i = 1:1:length(Location_disp)
    [x,y] = pol2cart(Location_disp(i,2),Location_disp(i,1));
    Location_cartX(i) = x;
    Location_cartY(i) = y;
    r = 0.001;
    % circle(x,y,r);
    % text(x,y,num2str(i));
    % grid on
end
Location_dispXY = [Location_cartX',Location_cartY'];
% axis equal
% title('Surface displacement and EM coil locations')


%% Roots of Bessel function of first kind, J_m(E) = 0, m= 0,1,2, ...

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


R = MdlParams.R;
lambda = E./R;  % Eqn 43, [Azhar 2008]

h = MdlParams.EMcoil2filmDist;
z = -h;
%% Compute Y and Z functions
Y = fnCompY(lambda, z, MdlParams); % function to compute Y, refer Eqn (35), [Azhar 2008]
Z = fnCompZ(lambda, z, MdlParams); % function to compute Z, refer Eqn (36), [Azhar 2008]
% Compute omega^2, F, Hmnc, Hmns, omega_d
omega2                              = fnCompOmega2(Y, Z, lambda, MdlParams); % function to compute omega^2, refer Eqn (50), [Azhar 2008]
Fmn                                 = fnCompF(lambda, Y, E, MdlParams); % function to compute Fmn, refer Eqn (52), [Azhar 2008]

%% Surface displacement as a function of coil input
for r_coils_in_count = 1:1:numel(r_coils_in_all)
    r_coils_in = r_coils_in_all(r_coils_in_count);
    theta_coils_in = theta_coils_in_all(r_coils_in_count);
    for r_in_count = 1:1:numel(r_in_all)
        r_in = r_in_all(r_in_count);
        theta_in = theta_in_all(r_in_count);

        
            % translation_coil_coord
            [Tx_coil,Ty_coil] = pol2cart(theta_coils_in,r_coils_in);
            new_coil_X = Tx_coil - Tx_coil;
            new_coil_Y = Ty_coil - Ty_coil;
            
            % translation_displacement_coord
            [Tx_disp,Ty_disp] = pol2cart(theta_in,r_in);
            new_in_X = Tx_disp - Tx_coil;
            new_in_Y = Ty_disp - Ty_coil;

            %calculate seperating distance
            dist = sqrt((new_coil_X - new_in_X)^2 + (new_coil_Y - new_in_Y)^2);


            [theta_coils_in_new,r_coils_in_new] = cart2pol(new_coil_X,new_coil_Y);
            [theta_in_new,r_in_new] = cart2pol(new_in_X,new_in_Y);
            

            % if dist <= Min_distance

                [H0_coils, Hmnc_coils, Hmns_coils]  = fnCompH_coils(lambda, E, r_coils_in_new, theta_coils_in_new); % function to compute H_mnc^j and H_mns^j,
                % refer Eqn (53)  and Eqn(56), [Azhar 2008]
                omega_d                             = fnCompOmegad(lambda, MdlParams); % function to compute omega_d_mn, refer Eqn (58), [Azhar 2008]
                [H0, Hmnc, Hmns]                    = fnCompH(lambda, r_in_new, theta_in_new); % function to compute H_mnc and H_mns, refer Eqn (44), [Azhar 2008]
                %% Compute A, B, C, D matrices

                A = fnCompA(omega2, omega_d); % function to compute system's A matrix, Pg 135 after Eqn(59), [Azhar 2008]
                B = fnCompB(Fmn, H0_coils, Hmnc_coils, Hmns_coils, length(r_coils_in_new)); % function to compute system's B matrix, Pg 135 & Pg 136after Eqn(59), [Azhar 2008]
                B_mimo(:,r_coils_in_count) =  B;
                C = fnCompC(H0, Hmnc, Hmns); % function to compute system's C matrix, Pg 136 after Eqn(59), [Azhar 2008]
                C_mimo(r_in_count,:) = C;
                D = zeros(length(r_in), length(r_coils_in_new)); % D matrix is 0 matrix
                %% Compute Transfer function, Bode plot and Step response
                N_coil = r_coils_in_count;
                N_disp = r_in_count;
                
                % transfer function
                G(r_in_count,r_coils_in_count) = tf(ss(A,B,C,D));                

    end
end

%% 
[Nrow, Ncol] = size(G);
G0           = zeros(Nrow, Ncol);
for ii=1:Nrow
    for jj=1:Ncol
        numG = cell2mat(G(ii,jj).Numerator);
        denG = cell2mat(G(ii,jj).Denominator);
        numGdc = numG(end);
        denGdc = denG(end);
        G0(ii, jj) = numGdc/denGdc;
    end
end
%% Mirror surface - Reference_signals

% figure(5);
x = linspace(min(Location_cartX.*1e3),max(Location_cartX.*1e3),201);
a = 34225;
y = x.^2/(4*a);
% plot(x,y)
% xlabel('mm')
% ylabel('mm')

r_in_all = r_in_all*1e3; %mm
for i = 1:1:length(r_in_all)
   Reference_signals(i,1) = interp1(x,y,r_in_all(i)) + 1e-4; % correction for porous media thickness
end


%% PI controller

kp = 2000; % 100
ki = 5; % 0.3

K_tf = tf([kp, ki],[1 0]);

I = eye(numel(r_coils_in_all));

G0inv = inv(G0);

K_bar = K_tf * eye(size(G0));

% Slopes and Start Times

db = 1;

m_tStart = [ db*1e-3 3;
             -db*1e-3 8.5];

m1      = m_tStart(1,1);
tStart1 = m_tStart(1,2);
m2      = m_tStart(2,1);
tStart2 = m_tStart(2,2);

%OL
m_tStart_OL = [ db*0.5e-3    3;
             -db*0.5e-3   8.5];

m1_OL      = m_tStart_OL(1,1);
tStart1_OL = m_tStart_OL(1,2);
m2_OL      = m_tStart_OL(2,1);
tStart2_OL = m_tStart_OL(2,2);

tSimStop = 8.5;

% disturbance gains
g1 = 0.0; g2 = 1; g3 = 0.5; g4 = 0.-0.5;  g5 = -1; g6 = -0.5; g7 = 0.5;
gain_vec_7pts = [g1,g2,g3,g4,g5,g6,g7]';


% Open loop input
u1 = 0.341009;
u27 = 0.209818;



%% Simulink
out = sim("MIMO_MFDM_Mdl_Zenith.slx");

%% Simulink post process - CL
out.logsout; % read log data from simulink

Output_Signals_CL = get(out.logsout,"Output_Signals_CL");
Output_Signals_data_7pts_CL = Output_Signals_CL.Values;

%% Simulink post process - OL

Output_Signals_OL = get(out.logsout,"Output_Signals_7pts_OL");
Output_Signals_data_7pts_OL = Output_Signals_OL.Values;



% control signal
Control_Signals_CL = get(out.logsout,"Control_Signals_CL");
Control_Signals_data_CL = Control_Signals_CL.Values;

figure(3)
plot(Output_Signals_data_7pts_CL.Time,Output_Signals_data_7pts_CL.Data(:,1)*1e3,'Linewidth',1.5)
hold on;
plot(Output_Signals_data_7pts_CL.Time,Output_Signals_data_7pts_CL.Data(:,2)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_CL.Time,Output_Signals_data_7pts_CL.Data(:,3)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_CL.Time,Output_Signals_data_7pts_CL.Data(:,4)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_CL.Time,Output_Signals_data_7pts_CL.Data(:,5)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_CL.Time,Output_Signals_data_7pts_CL.Data(:,6)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_CL.Time,Output_Signals_data_7pts_CL.Data(:,7)*1e3,'Linewidth',1.5)
grid on
title('CL Response of points above coils');
xlabel('Time (s)');
ylabel('Amplitude (um)');
xlim([0,8.5]);
legend('Point 1','Point 2','Point 3','Point 4','Point 5','Point 6','Point 7');

figure(4)
plot(Output_Signals_data_7pts_OL.Time,Output_Signals_data_7pts_OL.Data(:,1)*1e3,'Linewidth',1.5)
hold on;
plot(Output_Signals_data_7pts_OL.Time,Output_Signals_data_7pts_OL.Data(:,2)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_OL.Time,Output_Signals_data_7pts_OL.Data(:,3)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_OL.Time,Output_Signals_data_7pts_OL.Data(:,4)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_OL.Time,Output_Signals_data_7pts_OL.Data(:,5)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_OL.Time,Output_Signals_data_7pts_OL.Data(:,6)*1e3,'Linewidth',1.5)
plot(Output_Signals_data_7pts_OL.Time,Output_Signals_data_7pts_OL.Data(:,7)*1e3,'Linewidth',1.5)
grid on
title('OL Response of points above coils');
xlabel('Time (s)');
ylabel('Amplitude (um)');
xlim([0,8.5]);
legend('Point 1','Point 2','Point 3','Point 4','Point 5','Point 6','Point 7');


% figure(5)
% plot(Control_Signals_data_CL.Time, Control_Signals_data_CL.Data(:,1),'Linewidth',1.5)
% hold on;
% plot(Control_Signals_data_CL.Time, Control_Signals_data_CL.Data(:,2),'Linewidth',1.5)
% plot(Control_Signals_data_CL.Time, Control_Signals_data_CL.Data(:,3),'Linewidth',1.5)
% plot(Control_Signals_data_CL.Time, Control_Signals_data_CL.Data(:,4),'Linewidth',1.5)
% plot(Control_Signals_data_CL.Time, Control_Signals_data_CL.Data(:,5),'Linewidth',1.5)
% plot(Control_Signals_data_CL.Time, Control_Signals_data_CL.Data(:,6),'Linewidth',1.5)
% plot(Control_Signals_data_CL.Time, Control_Signals_data_CL.Data(:,7),'Linewidth',1.5)
% grid on
% title('Control signal of coil 1 and 2');
% xlabel('Time (s)');
% ylabel('Control currents (mA)');
% legend('Coil 1','Coil 2');


