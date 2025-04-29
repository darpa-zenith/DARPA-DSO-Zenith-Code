clc;clear all;close all;
load act35mm.mat
load xAct.mat
load yAct.mat
load currentdistribution.mat

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

xdata=xAct;
ydata=yAct;
nact=271;
ax=act35mm(:,1);
ay=act35mm(:,2);
actuator=act35mm;
np = size(xdata,2);
nloops = size(actuator,1);
actdia = 0.035; % actuator diameter
actoffset = -0.035; % actuator offset
xact = ax*actdia; % x position of the actuators on the grid
yact = ay*actdia; % y position of the actuators on the grid
zdata = 0.125*(xdata.^2 + ydata.^2); % parabolic surface of data
zact = 0.125*(xact.^2 + yact.^2); % offset of the actuators to maintain uniform distance from parabolic surface
infBx = zeros(np,nact); infBy = zeros(np,nact) ; infBz = zeros(np,nact);
EMcurrent = 0.1;
HelmholtzCurrent = 100;

% calculate magnetic field of Helmholtz coil
coil_top = [0.4 0.2]; % radius, height of 1st coil
coil_down = [0.4 -0.2]; % radius, height of 2nd coil
[bx_top,by_top,bz_top] = magloop(xdata,ydata,zdata,coil_top(1),coil_top(2));
[bx_down,by_down,bz_down] = magloop(xdata,ydata,zdata,coil_down(1),coil_down(2));
% add magnetic field components from 1st and 2nd coils
bx_coil = bx_top + bx_down;
by_coil = by_top + by_down;
bz_coil = bz_top + bz_down;

% compute influence matrix
for i = 1:nact % iterate for each actuator
    for j = 1:np % iterate through all data points
        for k = 1:nloops % iterate through all current loops
            zcoil = actuator(k,2) + zact(i) + actoffset;
            [dBx,dBy,dBz] = magloop(xdata(j)-xact(i),ydata(j)-yact(i),zdata(j),actuator(k,1),zcoil);
            infBx(j,i) = infBx(j,i)+dBx;
            infBy(j,i) = infBy(j,i)+dBy;
            infBz(j,i) = infBz(j,i)+dBz;
        end
    end
end

% physical constants
ur = 2.7;  % magnetic permeability of the ferrofluid         
rho = 1210;  % density of the ferrofluid      
u0 = 4*pi*1e-7; % magnetic permeability of vaccum
g = 9.8;      % acceleration gravity
k = (ur-1)/(u0*ur*rho*g);

% compute deformation
w = k*(ur*(infBx.*Current'+ bx_coil*HelmholtzCurrent).^2+ur*(infBy.*Current' + by_coil*HelmholtzCurrent).^2+(infBz.*Current' + bz_coil*HelmholtzCurrent).^2);

figure;
scatter3(xAct,yAct,w(1,:)*1e6,[],w(1,:)*1e6);
colormap jet
view(2)
xlabel('X coordinates (m)')
ylabel('Y coordinates (m)')
a=colorbar;
ylabel(a,'Deformation (um)','FontSize',16,'Rotation',270);