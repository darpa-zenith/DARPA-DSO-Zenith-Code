clc;clear all;close all;

% Function that computes the influence matrices of an array of coil actuators located inside an Helmholtz coil
%
% INPUTS:
%     actuator: Matrix containing the coil actuator geometry (nloops x 2)
%               where nloops is the number of loops that the coil is made of,
%               1st column is the ith loop radius, and 2nd column is the ith loop z position
%         nact: Number of actuators in the array (scalar)
%       ax, ay: x and y coordinates of the actuators on the array in unit of the actuator diameter (1 x nact vectors)
% xdata, ydata: x and y coordinates at which the influence is calculated (vectors)
%               * Use Matlab meshgrid function to create matrices then transpose to vectors
%               * Number of data points is np
% OUTPUTS:
%         infBx, infBy, infBz: Magnetic field influence functions (np x nact)  
%   bx_coil, by_coil, bz_coil: Values of the magnetic field components created by the Helmholtz coil (1 x np)  
%
% PREREQUISITES: magloop.m    
%
% initialiaze data

[X,Y] = meshgrid(linspace(-0.4,0.4,100));

xdata=X;
ydata=Y;
IdealCurrent=100;
% calculate magnetic field of Helmholtz coil
coil_top = [0.4 0.2]; % radius, height of 1st coil
coil_down = [0.4 -0.2]; % radius, height of 2nd coil
[bx_top,by_top,bz_top] = magloop(xdata,ydata,0,coil_top(1),coil_top(2));
[bx_down,by_down,bz_down] = magloop(xdata,ydata,0,coil_down(1),coil_down(2));
% add magnetic field components from 1st and 2nd coils
bx_coil = bx_top + bx_down;
by_coil = by_top + by_down;
bz_coil = bz_top + bz_down;

figure;
plot(X(50,:), IdealCurrent*bz_coil(50,:)*10000) % unit Gauss
xlabel('coordinates (m)')
ylabel('Magnetic Field (Gauss)')
