%% 2D Halbach Array Simulation
% Computes the magnetic flux density map, equipotential line contour map,
% as well as the kelvin body force distribution of a 2D flat halbach array
% as a function of the geometrical parameters. The code approximates the 
% magnets as infinite current sheets going into and out of the 
% side of the magnet.
% 
% 
% V1.2 - Hugh Chen (02/04/2024)
% V1.1 - Tianyang Hu (08/01/2024)
clc
clear
close all

%Path
addpath('./functions')

% Simulation parameters
w           = 1/2*0.0254;                      % [m]     Width of magnet
h           = 1/2*0.0254;                      % [m]     Height of magnet 
N_coil      = 20;                              % number of magnets
Ms          = 875;                             % Saturation Magnetization (A/m)
mu0         = 4*pi*1e-7;                       % [N/A^2] Permeability of 
% free space
Mm          = (1.48/mu0)*ones(N_coil,1);       % [A/m]   Magnetization
Ke          = Mm;                              % [A/m]   Sheet current         

% Simulation resolution parameter
gridsz      = 2000;                            % Grid size for contour plot
gridsz_q    = 50;                              % Grid size for quiver plot


%Storing things as param
params.w        = w;
params.h        = h;
params.N_coil   = N_coil;
params.Ms       = Ms;
params.mu0      = mu0;
params.Mm       = Mm;
params.Ke       = Ke;
params.gridsz   = gridsz;
params.gridsz_q = gridsz_q;

% magnet placement (x,y coordinates)
r0                          = zeros(N_coil,2);
r0(1:N_coil/2,1)            = w/2;
r0(N_coil/2+1:N_coil,1)     = -w/2;
r0(1:N_coil/2,1)            = r0(1:N_coil/2,1) - w*[N_coil/2:-1:1]';
r0(N_coil/2+1:N_coil,1)     = r0(N_coil/2+1:N_coil,1) + w*[1:1:N_coil/2]';

% magnet magnetization direction
theta           = zeros(N_coil,1);
theta(1:4:end)  = deg2rad(0);
theta(2:4:end)  = deg2rad(90);
theta(3:4:end)  = deg2rad(180);
theta(4:4:end)  = deg2rad(270);

% width input for calculating B field
w_vec           = ones(N_coil,1)*w;
w_vec(2:2:end)  = h;

% height input for calculating B field
h_vec           = ones(N_coil,1)*h;
h_vec(2:2:end)  = w;


rangex          = 1.8*max(r0(:,1));                   % X Range of the plot
rangey          = 0.8*max(r0(:,1));                   % Y Range of the plot
[X,Y]           = meshgrid(linspace(-rangex,rangex,gridsz),linspace(-rangey,rangey,gridsz)+0.01);
hx_eric         = 2*rangex/(gridsz - 1);              % resolution for gradient function in x dir.
hy_eric         = 2*rangey/(gridsz - 1);              % resolution for gradient function in y dir.
[Xq,Yq]         = meshgrid(linspace(-rangex,rangex,gridsz_q),linspace(-rangey,rangey,gridsz_q)+0.01);

%storing params
params.hx_eric  = hx_eric;
params.hy_eric  = hy_eric;
% Calculates B fields [T]

% B field for contour plot
[Bx,By]         = calculatingB(params,h_vec,w_vec,X,Y,theta,r0); 
% B field for quiver plot
[Bx_q,By_q]     = calculatingB(params,h_vec,w_vec,Xq,Yq,theta,r0);
% B field in magnitude
B               = sqrt(Bx.^2+By.^2);                                     
        

%% Kelvin Body Force Map

% Ferrofluid parameters, obtained rho and chi0 specs from Ms interpolation
[rho,chi0,phi] = get_Rho_Chi0_Phi(Ms); 
gamma          = 3*chi0/Ms;
%Storing params
params.rho     = rho;
params.chi0    = chi0;
params.phi     = phi;
params.gamma   = gamma;

Hx             = Bx / mu0;                       %obtaining X magnetic field
% from B field
Hy             = By / mu0;                       %obtaining Y magnetic field 
% from B field
H              = sqrt(Hx.^2+Hy.^2);
nx_H           = Hx./H;                          %unit vector of H in x dir
ny_H           = Hy./H;                          %unit vector of H in y dir
M              = magnetization(H,params);        %M as scalar
M_vec          = cat(3,M.*nx_H,M.*ny_H);         %M as vector form, as M is 
% in direction of H

[Fk_x,Fk_y,Fk] = Kelvin_H(M_vec, Hx,Hy,params); %Kelvin body force in component form
%and scalar, all in g

% quiver mesh
gridsz_Fk      = 50;
[Xq_Fk,Yq_Fk]  = meshgrid(linspace(-rangex,rangex,gridsz_Fk),linspace(-rangey,rangey,gridsz_Fk)+0.005);
Fk_x_q         = interp2(X,Y,Fk_x,Xq_Fk,Yq_Fk);     %interp between known kelvin forces in x dir
Fk_y_q         = interp2(X,Y,Fk_y,Xq_Fk,Yq_Fk);     %interp between known kelvin forces in y dir

%% Plot 2D Halbach Array with its B field
% Plot the B field strength contour
fig1           = figure;
contourf(X,Y,B,linspace(0,1,200))
c              = colorbar('FontSize',40);
c.Label.String = 'B field strength (Tesla)';
hold on

quiver(Xq,Yq,Bx_q,By_q,'Color',[1 1 1])
axis equal

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end

%% Plot B vector field
grid on
% title('B Field of Halbach Array')
xlabel('X (m)','FontSize',40)
ylabel('Y (m)','FontSize',40)


%% Plot 2D Halbach Array with equipotential lines 
% Equipotential
MagFun         = @(Ht) Ms*(coth(gamma*Ht)-1/gamma./Ht);
Pi_m           = @(Ht) -mu0/rho * trapz(1:Ht,MagFun(1:Ht));

% Plot the equipotential field contour
fig2           = figure;
contourf(X,Y,0.5*mu0*chi0*H.^2 / rho + 9.81*Y,linspace(0,1,200))
c              = colorbar('FontSize',40);
c.Label.String = 'Equipotential (m^2/s^2)';
hold on

quiver(Xq,Yq,Bx_q,By_q,'Color',[1 1 1])
axis equal

% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end


% Plot B vector field
grid on
% title('Equipotential lines')
xlabel('X (m)','FontSize',40)
ylabel('Y (m)','FontSize',40)




%% Plot the Kelvin Body Force contour map
% levels         = [linspace(0,120,20)];
levels         = [linspace(0,40,20)];
fig3           = figure;
contourf(X,Y,Fk,levels)
c              = colorbar('FontSize',40);
c.Label.String = 'Kelvin Body Force Strength (g)';
hold on

% Adding the quiver on the KBF map
quiver(Xq_Fk,Yq_Fk,Fk_x_q,Fk_y_q,'r')
axis equal
% Plot the magnetized material shapes
for i = 1:size(r0,1)
    Plot_magnets(theta(i),w_vec(i),h_vec(i),r0(i,:),1);
end

% title('Kelvin Body Force of Halbach Array')
xlabel('X (m)','FontSize',40)
ylabel('Y (m)','FontSize',40)

%% Imposing LGST lab plotting standard
set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
set(findobj('Type', 'axes'), 'FontName','latex','FontSize', 40,...
    'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);


%saving the plot in fullscale pdf
set(fig1, 'WindowState', 'maximized')
% exportgraphics(fig1, 'Magnetic_Flux_Density_2D_flat_HB.pdf', 'ContentType', 'vector');

exportgraphics(fig1, 'Magnetic_Flux_Density_2D_flat_HB.pdf', "Resolution",150);

set(fig2, 'WindowState', 'maximized')
% exportgraphics(fig2, 'EQP_line_map_2D_flat_HB.pdf', 'ContentType', 'vector');
exportgraphics(fig2, 'EQP_line_map_2D_flat_HB.pdf', "Resolution",150);

set(fig3, 'WindowState', 'maximized')
% exportgraphics(fig3, 'KBF_map_2D_flat_HB.pdf', 'ContentType', 'vector');
exportgraphics(fig3, 'KBF_map_2D_flat_HB.pdf', "Resolution",150);

