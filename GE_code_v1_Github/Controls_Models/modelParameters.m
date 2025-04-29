%% All model parameters are listed in this script

%%%%%% These parameters are obtained from Dr. Azhar Iqbal's thesis "Modeling and control 
%%%%%% of a Magnetic Fluid Deformable Mirror for Opthalmic Adaptive Optics Systems

% Conversion Factors
cm2m = 1e-2;
mm2m = 1e-3;
deg2rad = pi/180;


% Physical Constants
MdlParams.g                 = 9.81; % gravity [m/s^2]
%% Azhar Model Parameters
MdlParams.R                 = 100 * mm2m; % radius of liquid film layer [m]
MdlParams.thickness         = 2 * mm2m; % thickness of liquid film [m]
MdlParams.EMcoil2filmDist   = 1 * mm2m; % space between point object coil and bottom of liquid film [m]
MdlParams.HHcoilB0          = 17e-3; % Helmholtz coil uniform magnetic field strength [Tesla] or [N m^-1 A^-1], where A is ampere
MdlParams.mu0               = 4*pi*1e-7; % Magnetic permeability of Free space [H m^-1] or [N A^-2], where H is Henry and A is ampere
% Properties of ferrofluid EFH1
% Source: https://ferrofluid.ferrotec.com/products/ferrofluid-educational-fluid/efh/efh1/
% eta_factor                  = 25;%
eta_factor                  = 200;%
MdlParams.eta               = 5.8e-3 * eta_factor; % Dynamic viscosity [Pa-s] or [N-s/m^2] 

MdlParams.rho               = 1.21e3; % Density [kg/m^3] 
MdlParams.rel_mu            = 1.78; % Relative magnetic permeability [unitless]
MdlParams.mu                = MdlParams.rel_mu * MdlParams.mu0; % Magnetic permeability of fluid [H m^-1]
MdlParams.chi               = MdlParams.rel_mu - 1; % Magnetic susceptibility [unitless]
MdlParams.sigma             = 29e-3;% * sigma_factor; % Coefficient of surface tension [N/m] 