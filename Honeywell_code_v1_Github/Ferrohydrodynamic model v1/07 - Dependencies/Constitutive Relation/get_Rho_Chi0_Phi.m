function [rho,chi0,phi] = get_Rho_Chi0_Phi(Ms)
% interpoation via EMG water-based ferrofluid
% Calculate the diluted ferrofluid Initial Susceptibility chi0, 
% Magnetic Particle concentration phi, and Density rho [kg/m^3]
% from Saturation Magnetization Ms [A/m], based on interpolation
% of the EMG series water-based 
% ferrofluid physical property data.
%
% [rho,chi0,phi] = get_Rho_Chi0_Phi(Ms)
%
% Data for fit:
%     1st Input  : Ms
%
% Output:
%     1st Output : rho
%     2nd Output : chi0
%     3rd Output : phi
% V1.1 - Tianyang Hu(08/01/2024)


%% Interpolation data for EMG water-based series ferrofluid (Ferrotec Inc.)
% URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/
% EMG  0, 1:10-700, 1:5-700, 304, 308 408 507  508  509 605  607  700   705  707  708
% Initial Susceptibility (SI, vol)
Chi0_data = [0    0.181  0.39   5.03  0.5  0.5  1.63 0.88 0.5  3.02 1.63 12.57 4.04 1.51 0.63];               
% Saturation Magnetization [A/m]
Ms_data   = [0    3.971  5.2276 27.5  6.6  6.6  11   6.6  3.3  22   11   35.5  22   11   6.6]* 795.7747;      
% Magnetic Particle concentration [%]
phi_data  = [0    0.58   1.16   4.5   1.2  1.2  2    1.2  0.6  3.9  2    5.8   3.9  2    1.2];                
% Density @ 25 C [kg/m^3]
rho_data  = [1000 1020   1058   1240  1060 1070 1120 1070 1030 1180 1100 1290  1190 1100 1080];               


% Interpolators
[Phi_Ms_cfit, ~]   = createFit_PhiVsMs(Ms_data, phi_data);
[Chi0_phi_cfit, ~] = createFit_Chi0VsPhi(phi_data, Chi0_data);
[rho_phi_cfit, ~]  = createFit_RhoVsPhi(phi_data, rho_data);

% Values from interpolation
phi     = Phi_Ms_cfit(Ms);
rho     = rho_phi_cfit(phi);               % Density @ 25 C [kg/m^3]       
chi0    = Chi0_phi_cfit(phi);              % Saturation Magnetization [A/m]
end 

