%% Sample script for Constitutive relation equations
% This code demonstrates the usage of the constitutive relation functions
% to generate the a series of ferrofluid properties
%
% V1.0 - Hugh Chen (05/04/2024)

clear all

Ms             = 4000 %Example saturation magnetization, [A/m]
[rho,chi0,phi] = get_Rho_Chi0_Phi(Ms);