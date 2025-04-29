%% Passive and Active COMOSL simulation: Generate passive corrector shape for Comsol model
%
% The simulation procedure integrates both passive and active correctors 
% through a coupled MATLAB and COMSOL approach. Initially, the geometry 
% of the passive corrector, along with parameters such as the Halbach array
% dimensions, is defined in the MATLAB code. These parameters are then 
% used to update the COMSOL model, which begins the simulation process.
%
% Upon completion of the simulation, COMSOL sends the resulting data back 
% to MATLAB for further processing. In MATLAB, the overall potential is 
% evaluated to generate the EQP interface. Additionally, a Fast Fourier 
% Transform (FFT) analysis is performed on this data. 
%
% V1.0 Hugh (5/28/2024)


close all
plotting = true;
addpath('./functions')

% Inputs 
a0           = 0.000500238418579102; % Constant passive corrector thickness (m)
mu_passive   = 4.01558649902344; % Permeability of the passive corrector material
mu_active    = 3;                % Permeability of the active corrector material


%% Additional Geometry
h            = 9e-3;             % Ferrofluid height (over the halbach array)
delta        = 2e-3;             % Thickness of the ferrofluid layer
ME_thickness = 1e-3;             % Thickness of the ME blanket layer
w            = 0.5*0.0254;     % Width of the magnet (m)
b            = 1*0.0254;       % Depth of the magnet (m)
alpha        = 0;              % Tilting (rad)
xf           = 0.25;           % End of fluid interface (m)
ME_length    = 0.7;            % Length of magnetoelectric blanket (m)

% Ferrofluid parameters
fluct        = 0;                % Ferrofluid density fluctuation (0 to 1)
Ms           = 5607;             % Saturation Magnetization (of the Ferrofluid)
[rho,chi0,~]     = get_Rho_Chi0_Phi(Ms);  % Interpolated Ferrofluid data
FF_fluct_density = fluct;                 % fluctation in FF density
% Linear fits of the langevin curve of MH, empirical [Needed because of the
% particle migration problem]
M_linear_slope   = 9.286420326826829e-4;
M_linear_b	     = 3.929750472154023e3; 

% Constants
mu0              = 4*pi*1e-7;      % Magnetic permeability of vacuum (N/A^2)
g                = 9.81;           % Gravity (m/s2)

% Comsol file path
cs_file          ='./Halbach_WFE_200_um_improved.mph';
name             = 'GlobalWFE_active_with_passive_corrector';              %name for the save .m param

% Derived parameters
lambda           = 4*w;            % Magnetic wavelength (m)
k                = 2*pi/lambda;    % Wave number (rad/m)
A0               = 2*xf*delta;     % Interface area (m^2)

% Storing parameters
params.mu0       = mu0;
params.chi0      = chi0;
params.Ms        = Ms;
params.rho       = rho;
params.w         = w;
params.delta     = delta;
params.b         = b;
params.h         = h;
params.lambda    = lambda;
params.k         = k;
params.g         = g;
params.alpha     = alpha;
params.xf        = xf;
params.A0        = A0;
params.a         = a0;
params.FF_fluct_density = FF_fluct_density;
params.ME_length        = ME_length;
params.mu_passive       = mu_passive;
params.ME_thickness     = ME_thickness;
params.mu_active        = mu_active;

%% Magnetic passive corrector
% Discretization
x_ms           = [-0.35:2e-4:0.35];
% surface profile:
y_ms           = a0 .* ones(1,length(x_ms)) - h + ME_thickness;
% safety conditions:
y_ms(y_ms<-h)  = -h + ME_thickness; 
y_ms(1)        = -h + ME_thickness;
y_ms(end)      = -h + ME_thickness; 

% Storing parameters
params.x_ms    = x_ms;
params.y_ms    = y_ms;

%% Run COMSOL model
if ~exist('model')
    model = mphopen(cs_file);
end

% Update parameters for COMSOL model
model.param().set('w',[num2str(w),' [m]']);                         % Updating the magnet width in COMSOL
model.param().set('b',[num2str(b),' [m]']);                         % Updating the magnet height in COMSOL
model.param().set('hff',[num2str(h),' [m]']);                       % Updating ferrofluid height
model.param().set('delta',[num2str(delta),' [m]']);                 % Ferrofluid thickness, one potential new variable
model.param().set('Ms',[num2str(Ms),' [A/m]']);                     % Saturation magnetization of the ferrofluid
model.param().set('FF_fluct_density',[num2str(FF_fluct_density)]);  % Ferrofluid density fluctuation
model.param().set('k',[num2str(k),' [1/m]']);
model.param().set('ME_thickness',[num2str(ME_thickness)]);
model.param().set('ME_length',[num2str(ME_length)]);
model.param().set('mu_passive',[num2str(mu_passive)]);
model.param().set('mu_active',[num2str(mu_active)]);
model.param().set('M_linear_slope',[num2str(M_linear_slope)]);
model.param().set('M_linear_b',[num2str(M_linear_b),' [A/m]']);

% Update geometry
model.geom('geom1').feature('pol1').set('table', [x_ms',y_ms']); % Updating the passive corrector geometry in COMSOL

% Rebuild geometry [important to ensure materials and BC are correctly
% assigned]
model.component('comp1').geom('geom1').run;

% Reset material and physics
model.material('mat7').selection.set([8]); % be aware of selected component

% Run study
model.study('std1').run;


% Extract data
% limiting the x range to avoid meshing prob near FF bound:
x_ff        = (-xf+0.01):1e-3:(xf-0.01);

% ensuring that the y_ff is having a sufficiently high fidelity while
% matching vectors ensures residuals_setPI work:
y_ff        = linspace(delta/4,delta-delta/4,length(x_ff));
[X_ff,Y_ff] = ndgrid(x_ff,y_ff);

H  = mphinterp(model,'mf.normH','coord',[X_ff(:),Y_ff(:)]'); %using COMSOL data to interpolate normH data
Mx = mphinterp(model,'mf.Mx','coord',   [X_ff(:),Y_ff(:)]'); %using COMSOL data to interpolate Mx data
My = mphinterp(model,'mf.My','coord',   [X_ff(:),Y_ff(:)]'); %using COMSOL data to interpolate My data

% Convert data (rearrange them accordingly)
[~,~,H]     = spacingconv(X_ff,Y_ff,X_ff(:),Y_ff(:),H);
[~,~,Mx]    = spacingconv(X_ff,Y_ff,X_ff(:),Y_ff(:),Mx);
[~,~,My]    = spacingconv(X_ff,Y_ff,X_ff(:),Y_ff(:),My);


% Evaluating Magnetic potential (m2/s2)
Pi_m        = gridinterp_Pi_m(X_ff,Y_ff,H,params);
% Evaluating Gravity potential (m2/s2)
Pi_g        = griddedInterpolant(X_ff,Y_ff,g*Y_ff,'linear');

% Storing parameters
params.x_ff = x_ff;
params.y_ff = y_ff;
params.X_ff = X_ff;
params.Y_ff = Y_ff;
params.H    = H;
params.Mx   = Mx;
params.My   = My;
params.Pi_m = Pi_m;
params.Pi   = @(x,y) (Pi_m(x,y) + Pi_g(x,y));

%% Determine equipotential line

% The target deviation provided by Optical team
p_target          = [-7127.63607978300	5293.05578812901	6500.05673636707	-2904.55896817343	-240.969106733449	27.0176476468912];
target_prof_x     = linspace(-0.25,0.25,length(x_ff));
target_prof_y     = polyval(p_target,target_prof_x)/1e6;

% Solving for the location of EQP line
options           = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-12,'OptimalityTolerance',1e-10,'StepTolerance',1e-10);

% Initial condition, just start from y = 0:
Y_IC              = ones(1,length(y_ff))*delta/2;
% Finding the f that minimize the difference in potential:
Y_t               = fsolve(@(Y)residuals_setPI(Y,params), Y_IC, options);

Y_to              = Y_t;
Deviation         = Y_to - Y_IC; % The distance from the ideal position


% Constructing the perfect passive corrector to serve as compare reference
p2                = polyfit(x_ff,Deviation,5);
passive_damping   = 0.01;
New_new_deviation = polyval(p2,x_ff) + passive_damping*(Deviation -polyval(p2,x_ff)) ;

%% EQP interface plotting
fig2      = figure;
if plotting
    hold on
    plot(target_prof_x,target_prof_y)
    plot(x_ff, Deviation)
    % plot(x_ff, New_new_deviation)
    xlabel('Location along FF layer [m]')
    ylabel('Potential Height [m]')
    legend('target','Passive and Active control')
end

% Maximum absolute discrepancy between real passive corrector + 100x passive corrector:
max_diff                 = max(abs(Deviation-New_new_deviation));

if max_diff == 0
    max_diff = 1;
end

params.max_diff          = max_diff;
params.Y_to              = Y_to;
params.target_prof_x     = target_prof_x;
params.target_prof_y     = target_prof_y;
params.Deviation         = Deviation;
params.New_new_deviation = New_new_deviation;

%% FFT plotting
%FTT part
n           = length(x_ff);                % the length of the eqp vector
% Single-Sided Amplitude Spectrum of equipotential line
fs          = 1/(x_ff(2)-x_ff(1));         % the step size of the x dir
Y_tfft      = fft(Deviation - mean(Deviation));        % Computing the FFT of the signal (doing Y_t-mean(Y_t) to extract the deviation)
P2          = abs(Y_tfft/n);               % Rescaling the FFT result by n (since FFT has built in scale up) plus take complex magnitude of FFT
% Convert to single sided spectrum
P1          = P2(1:n/2+1);                 % Take the first half of the two-sided spectrum P2
P1(2:end-1) = 2*P1(2:end-1);               % Multiply the spectrum in the positive frequencies by 2, no need P1(1) and P1(end) by 2 because these amplitudes = zero & Nyquist frequencies
f           = fs/n*(0:(n/2));              % Define the frequency domain f for the single-sided spectrum

% Now plot the centered at 0 PSD
fig3        = figure(3);
hold on
plot(2*pi*f/k, 1e6*P1)
xlim([0,5])
ylim([0,100])
xlabel('Spatial Frequency [k]')
ylabel('Potential Amplitude [$\mu$m]')


%% Implementing Generic Format for plotting
if plotting
    set(findobj('Type', 'Legend'), 'Interpreter', 'latex', 'box','off');
    set(findobj('Type', 'axes'), 'FontName','Arial','FontSize', 16,...
        'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex','TickDir','out');
    set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
    set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [0    0.7620   16   11]); % [0 0 16 12]
    set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
    set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);
end

%% KBF calculation in the FF region (matching 90% of no-passive corrector)
Fk          = Kelvin_B(X_ff,Y_ff,params);
params.Fk   = Fk;
avgFk       = mean(mean(Fk));
params.avgF = avgFk;


%% Save .svg, and data
saveas(fig2,[name,'.svg'])
saveas(fig3,[name 'FFT','.svg'])
save([name,'-data.mat'],'params');
