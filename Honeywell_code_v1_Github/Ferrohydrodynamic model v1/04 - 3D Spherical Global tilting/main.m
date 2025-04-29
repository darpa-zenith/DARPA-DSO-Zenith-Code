%% Analytical Halbach Tilting approximation
% Computes the ferrofluid interface as a function of the tilting angle and
% the liquid parameters. The code employs the surface rougness estimation
% obtained from the FHD Bernoulli equation and the first two modes of the
% Mallinson/Halbach discretization. This is superposed to the equipotential
% line with best spherical fit under tilting. Key missing elements: (i)
% magnetic edge effects, (ii) surface tension edge effects, (iii)
% cross-terms and second-order effects. The approximation is not valid for
% h/lambda < 0.1. 
%
% REFERENCE SYSTEM:
%
%  + Origin fluid circle (x = 0, z = 0)
%  |
%  |
%  | R       + Origin Halbach Array circle (x = 0, z = dO)
%  |         |
%  |         |
% \|/        |
%--'         |
%    delta   | Rm
%---         |
%    h       |
%---         |
%    b      \|/
%---         '
%
% V1.0 - Álvaro Romero-Calvo (05/04/2024)

clear all

% Path
addpath('./functions')

% Fluid parameters
delta  = 0.002;            % Fluid height (m)
sigma  = 61.7 * 1e-3;      % Surface tension (N/m)
alpha  = 10 * pi/180;      % Plate tilting (rad)
beta   = 10 * pi/180;      % Plate tiping (rad)
Ms     = 5607;             % Saturation magnetization ferrofluid (A/m)
R      = 2;                % Spherical mirror surface curvature radius (m)
A      = 0.5;              % Spherical mirror surface aperture (m)
DNoN   = 0;              % Percentual variation in ferrofluid density (-)

% Halbach parameters
w      = 1/2*0.0254;        % Width of individual magnets (m)
t      = 0;                 % Separation between magnets (m)
b      = 1*0.0254;          % Height of individual magnets (m)
L      = 0.30;              % Half-length of Halbach array (m)
h      = 9e-3;              % Separation between Halbach and FFL (m)
M0     = 1177.7e3;          % Magnet magnetization (A/m)
dO     = -0.077;               % Vertical displacement of the center of the array (m) - to compensate for gravity
Rm     = R + delta+h+b/2 + dO; % Magnet array radius at the centerline of the array (m)
Pc     = 1e-4;                 % Passive corrector (-);

% Magnetic array derived parameters
lambda = 4 * (w + t); % m
k      = 2*pi/lambda;

% Constants
mu0    = 4*pi*1e-7;           % Magnetic permeability of vacuum (N/A^2)
g      = 9.81;                % Gravity (m/s2)

% Computation
N      = 50;                  % Discretization points (2x in x, 1x in y, 4x in z)

% Constitutive equation
[rho,chi0,phi] = get_Rho_Chi0_Phi(Ms);   % Obtain rho (kg/m3), chi0, and phi0 (%vol) specs from Ms interpolation

% Store parameters
params.delta = delta;
params.rho   = rho;
params.phi   = phi;
params.alpha = alpha;
params.beta  = beta;  
params.sigma = sigma;
params.R     = R;
params.Rm    = Rm;
params.dO    = dO;
params.chi0  = chi0;
params.Ms    = Ms;
params.w     = w; 
params.t     = t;
params.b     = b;
params.L     = L;
params.M0    = M0;
params.g     = g;
params.mu0   = mu0;
params.h     = h;
params.lambda= lambda;
params.k     = k;
params.PC    = Pc;
params.A     = A;
params.DNoN  = DNoN;

%% Magnetic & reference computations
% Mesh (origin in focal point of R radius optical sphere)
[x_m, y_m, z_m] = ndgrid(linspace(-1.1*L,1.1*L,2*N),linspace(-1.1*L,1.1*L,2*N),linspace(-R - delta,-1.95,2*N));

% Evaluation grid for ferrofluid surface
x_t          = x_m(:,:,1);
y_t          = y_m(:,:,1);


% Dish profile
z_dish       = -sqrt((R + delta)^2- x_t.^2 - y_t.^2);

% Target ferrofluid profile (removing fluid outside of aperture area)
z_tar        = -sqrt(R^2 - x_t.^2 - y_t.^2);
outregion    = sqrt(x_t.^2+y_t.^2) > A/2;
z_tar(outregion) = z_dish(outregion);

% Ferrofluid half volume (to be conserved)
V0           = quad2d(@(x,y) interp2(x_t',y_t',z_tar'-z_dish',x,y,'linear'),...
               min(min(x_t)),max(max(x_t)),min(min(y_t)),max(max(y_t)),'AbsTol',1e-8);

% Coordinate transformation
z_p          = (Rm-b/2-h) - sqrt(x_m.^2 + y_m.^2 + (z_m - dO).^2);
x_p          = (Rm-b/2-h) * atan(-x_m./(z_m - dO));
y_p          = (Rm-b/2-h) * atan(-y_m./(z_m - dO));

% Definition of magnetic field (A/m)
H0           = M0 * exp(-k*h) * (1-exp(-k*b));
Hx           = H0 * exp(-k*z_p) .* sin(k*x_p);
Hy           = 0;
Hz           = H0 * exp(-k*z_p) .* cos(k*x_p);
H            = sqrt(Hx.^2 + Hy.^2 + Hz.^2);

% Magnetic potential (m2/s2)
Pi_m         = gridinterp_Pi_m(x_m,y_m,z_m,H,params);

% Gravity potential (m2/s2)
Pi_g         = griddedInterpolant(x_m,y_m,z_m,g*z_m,'linear');

% Save magnetic params
params.x_m  = x_m;
params.y_m  = y_m;
params.z_m  = z_m;
params.x_t  = x_t;
params.y_t  = y_t;
params.z_tar= z_tar;
params.z_dish= z_dish;
params.V0   = V0;
params.x_p  = x_p;
params.y_p  = y_p;
params.z_p  = z_p;
params.Pi_m = Pi_m;
params.Pi_g = Pi_g;
params.Pi   = @(x,y,z) (Pi_m(x,y,z) + Pi_g(x,y,z));

%% Calculate surface differential
% The Halbach array radius Rm is computed as the approximate solution
% imposed to the array in order to produce a spherical ferrofluid surface 
% with curvature R. Because such ideal profile is actually non-spherical
% and this first-order code is defined for a spherical substrate, the
% difference is subtracted from all solutions to facilitate the derivation
% of active control strategies. In this section, the equilibrium interface
% is computed at Zenith and the target profile is subtracted, providing the
% residual difference. 

% Integration options
opt       = optimoptions('fsolve','Display','iter');

% Zenith interface
z_z       = fsolve(@(Z)residuals_setPI(Z, -R, params), -sqrt(R^2 - x_t.^2 - y_t.^2), opt);

% Remove outer region
outregion    = sqrt(x_t.^2+y_t.^2) > A/2;
z_z(outregion) = z_dish(outregion);

% Differential height
dz        = z_z - z_tar;

% Save for upcoming correction
params.dz = dz;

%% Interface profile under tilting
% Update gravity potential under tilting
Pi_g        = griddedInterpolant(x_m,y_m,z_m,g*(- sin(alpha)*cos(beta)*x_m - sin(beta)*y_m + cos(alpha)*cos(beta)*z_m),'linear');
params.Pi_g = Pi_g;
params.Pi   = @(x,y,z) (Pi_m(x,y,z) + Pi_g(x,y,z));

% Integration options
opt       = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-10,'StepTolerance',1e-10);

% Interface shape - surface differential
z0          = fsolve(@(z0)residuals_setV(z0, params), -R, opt);
z_t         = fsolve(@(Z)residuals_setPI(Z, z0, params), -sqrt(R^2 - x_t.^2 - y_t.^2), opt) - dz;
 
% Remove outer region
z_t(outregion) = z_dish(outregion);

%% Surface mid-freq oscillations
dR           = R + delta - sqrt(x_t.^2 + y_t.^2 + z_t.^2); % Liquid depth (base of the plate to the FF surface)
epsilon      = Pc * midfreq_oscillations(x_t,abs(dR),params);
z_t          = z_t + epsilon .* cos(2*k*x_t);

% Clean external boundary
z_t(sqrt(x_t.^2+y_t.^2) > A/2) = NaN;

%% PLOTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represent interface and substrate
figure, hold on

% Interface
surf(x_t, y_t,z_t)

% Dish
z_dish2 = -sqrt((R + delta)^2- x_t.^2 - y_t.^2);
z_dish2(sqrt(x_t.^2+y_t.^2) > A/2+0.02) = NaN;
mesh(x_t, y_t,z_dish2,'EdgeColor','r') 

% Halbach array
z_hb = -sqrt((Rm-b/2)^2- x_t.^2 - y_t.^2)+dO;
z_hb(sqrt(x_t.^2+y_t.^2) > L) = NaN;
mesh(x_t, y_t,z_hb,'EdgeColor','k') 

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
cb = colorbar;
ylabel(cb,'z(m)')
axis equal
legend('Ferrofluid interface', 'Dish', 'Halbach array surface','location','southeast')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represent interface and substrate in 2D
figure, hold on

% Parameters
line = N;

% Interface
plot(x_t(:,line),z_t(:,line),'b')

% Dish
plot(x_t(:,line),z_dish2(:,line),'r')

% Halbach array
plot(x_t(:,line),z_hb(:,line),'k')

xlabel('x (m)')
ylabel('z (m)')
axis equal % change to axis normal for an expanded view
legend('Ferrofluid interface', 'Dish', 'Halbach array surface')

%% Write data to text file
writematrix(x_t,sprintf('./data/Interface_Ms%dA_m_lambda%din_h%dmm_%ddeg_x.txt',[Ms,lambda/0.0254, h*1000, alpha*180/pi]))
writematrix(y_t,sprintf('./data/Interface_Ms%dA_m_lambda%din_h%dmm_%ddeg_y.txt',[Ms,lambda/0.0254, h*1000, alpha*180/pi]))
writematrix(z_t,sprintf('./data/Interface_Ms%dA_m_lambda%din_h%dmm_%ddeg_z.txt',[Ms,lambda/0.0254, h*1000, alpha*180/pi]))

%% Plotting standardization
set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
    'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [0    0.7620   16   11]); % [0 0 16 12]
set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);