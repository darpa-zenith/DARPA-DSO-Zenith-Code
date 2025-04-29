%% 3D Halbach Array COMSOL Simulation
% The script first load a pre-built COMSOL model "HB_Environment", then it
% will add geometry, materials, physics, mesh to build a 3D simulation of
% a parabolic Halbach Array with/without a magnetic Passive Corrector (Soft
% Iron without losses) above it. The simulation include the magnetic flux density field 
% (B field) created by the Halbach Array, ferrofluid body force potential field, 
% electromagnetic force and torque experienced by each BMN-52 magnet. 
% Finally, it will extract the ferrofluid body force equipotential lines 
% from the 6 cutplanes in the model and plot them with the ideal surface 
% shape. It will also extract force data from the COMSOL model and generate
% a 2D force map. 
%
% V 1.0, Tianyang Hu, 6/23/2024
% Email: thu98@gatech.edu

%%
clc;clear;close all
addpath("functions/","data/")
model = mphload('HB_Environment.mph');                                            % Load the pre-build environment model

%% Initialize the design parameters [User Defined]
tic
% Define the Square Halbach array geometry (lengths are equal in x, and y dimensions )
w                  = 1 * 0.0254;                                                  % magnet width (m) 
d                  = w/20;                                                        % gap between magnets (m) 
n_mag              = 28;                                                          % Number of magnets in one dimension (assumed to be multiple of 4)
half_range         = (n_mag-1)/2 * (d+w);                                         % The last magnets location in half range
% Define the parabolic shape of the Halbach Array
a                  = 0.266;                                                       % x dimension curvature z = a*x^2
b                  = 0.266;                                                       % y dimension curvature z = b*y^2
model.param.set('a', num2str(a));
model.param.set('b', num2str(b));

% Define a Halbach Array to ferrofluid height (Assuming Halbach Array passes the origin [0, 0, 0])  
h_ff2hb            = 0.009;                                                       % Ferrofluid bottom to Halbach Array top distance[m]
h_pc2hb            = 0.003;                                                       % Passive Corrector bottom to Halbach Array top distance [m] 
delta_ff           = 0.002;                                                       % Ferrofluid thickness [m]


% Define a 3D cutpoint position where the equipotential line will be
% evaluated using the potential value at this cutpoint, default position is
% the x = 0, y = 0, z = ferrofluid top surface
cutpoint           = [0;0;h_ff2hb+w/2+delta_ff];

% Define passive corrector geometry
A                  = 0.002;                                                        % Amplitude   [m]
eps                = 0.002;                                                        % Minimum thickness   [m] 
lambda_hb          = 2*(w+d);                                                      % Halbach array Wavelength  [m]
lambda_pc_bottom   = lambda_hb-2*sinh(n_mag*lambda_hb*b)*h_pc2hb/n_mag;            % Passive Corrector Wavelength [m]
L_hb               = n_mag/2*lambda_hb;                                            % Total Halbach Array length [m]

% Switches of Passive Corrector and force calculation
use_pc             = 0;                                                            % 1 for adding a Passive Corrector
calculate_force    = 1;                                                            % 1 for adding force calculation nodes for magnets

% Choose 'Step' for only moving magnets up and down, 'Rotate' for moving
% and rotating magnets such that the top surface of magnets aligned smoothly. 'Rotate' option requires more memory
StepOrRotate       = 'Rotate';                                  
 switch StepOrRotate
     case 'Step'
         [r_mag,Br_orient_strArr] = StepHB_magnetPose(a,b,n_mag,half_range);                % Calculate the magnets positions and magnetization orientation (Step Halbach Array)
     case 'Rotate'
         [r_mag,k,psi,Br_orient_strArr] = RotatingHB_magnetPose(a,b,n_mag,half_range);      % Calculate the magnets positions, single-rotation axis (cartesian), single-rotation angle (deg) and magnetization orientation (Rotating  Halbach Array)
 end

%% Build the Geometries
% Set the dimension of the air block centered at [0,0,0] [m] (The simulation boundary should be large enough to be approximated as inifinity)
model.geom('geom1').feature('blk1').set('size', [6;6;6]);           % [m]

% Build Halbach Array geometry
f = waitbar(0, 'Starting');                                         % Initialize the waitbar
for i = 1:n_mag*n_mag                                               % For each magnet
    if i == 1                                                       % Create 4 Cumulative Selection nodes, each cumulative selection will contain the magents whose magnetization orientations are the same (up, left, down, right)
        model.component('comp1').geom('geom1').selection.create('csel1', 'CumulativeSelection');
        model.component('comp1').geom('geom1').selection.create('csel2', 'CumulativeSelection');
        model.component('comp1').geom('geom1').selection.create('csel3', 'CumulativeSelection');
        model.component('comp1').geom('geom1').selection.create('csel4', 'CumulativeSelection');
    end
    % Add geometry nodes
    model.geom('geom1').create(['blk',num2str(i+1)], 'Block');                                                                   % Create a geometry node
    model.geom('geom1').feature(['blk',num2str(i+1)]).set('size', [w;w;w]);                                                      % Define the magnet dimension
    model.geom('geom1').feature(['blk',num2str(i+1)]).set('pos', r_mag(i,:)');                                                   % Define the magnet positon 
    model.geom('geom1').feature(['blk',num2str(i+1)]).set('base', 'center');                                                     % Define magnet centers as the position
    model.geom('geom1').feature(['blk',num2str(i+1)]).set('selresult', true);                                                    % Enable the selection  
    model.geom('geom1').feature(['blk',num2str(i+1)]).set('contributeto', ['csel',num2str(mod(i-1,4)+1)]);                       % Assign geometry into one of the four Cumulative Selection nodes
    % Add Fillet nodes
    model.component('comp1').geom('geom1').create(['fil',num2str(i+1)], 'Fillet3D');                                                          % Create a fillet node
    model.component('comp1').geom('geom1').feature(['fil',num2str(i+1)]).selection('edge').set(['blk',num2str(i+1)], [1 2 3 6 7 9 10 12]);    % Select the 12 edges of the cube magnet in fillet node
    model.component('comp1').geom('geom1').feature(['fil',num2str(i+1)]).set('radius', '0.0005');                                             % Assign a radius of curvature to the fillet.
    
    if strcmp(StepOrRotate,'Rotate')
        % Add Rotation nodes
        model.component('comp1').geom('geom1').create(['rot',num2str(i+1)], 'Rotate');                                                        % Create a rotation node
        model.component('comp1').geom('geom1').feature(['rot',num2str(i+1)]).selection('input').set({['fil',num2str(i+1)]});                  % Select the geometry in rotation node 
        model.component('comp1').geom('geom1').feature(['rot',num2str(i+1)]).set('axistype', 'cartesian');                                    % Set the rotation axis type as a cartesian vector
        model.component('comp1').geom('geom1').feature(['rot',num2str(i+1)]).set('ax3', k(i,:));                                              % Set the single-rotation axis
        model.component('comp1').geom('geom1').feature(['rot',num2str(i+1)]).set('rot', psi(i));                                              % Set the single-rotation angle
        model.component('comp1').geom('geom1').feature(['rot',num2str(i+1)]).set('pos', r_mag(i,:));                                          % Set the position of the origin of rotation axis (need to be the same as the magnets position) 
    end

    % Update waitbar
    waitbar(i/(n_mag*n_mag), f , sprintf('Geometry Progress: %d %%', floor(i/(n_mag*n_mag)*100)));
end

% Build geometries
model.geom('geom1').run

% Build Magnetic Passive Corrector geometry
if use_pc
        model.param.set('A', num2str(A));
        model.param.set('lambda', num2str(lambda_pc_bottom));
        model.param.set('eps_shim', num2str(eps));
         
        d_bottom = h_pc2hb+w/2;
        model.param.set('d_bottom', num2str(d_bottom));
        d_top        = '(A*abs(sin(((2*b*y_top*sqrt(1+4*b^2*y_top^2)+asinh(2*b*y_top))/(4*b))*2*pi/lambda))*((mod(y_top,lambda)>lambda/2) && (mod(y_top,lambda)<lambda))+eps_shim)';
        
        x_pc_bottom_min = -(half_range+w/2);                                                                
        x_pc_bottom_max = half_range+w/2;           
        y_pc_bottom_min = -(half_range+w/2);
        y_pc_bottom_max = half_range+w/2;
        x_pc_top_min = -(half_range+w/2);
        x_pc_top_max = half_range+w/2;
        y_pc_top_min = -(half_range+w/2);
        y_pc_top_max = half_range+w/2;
             
        % Define the expression for the bottom, top surfaces of the Passive Corrector
        z_bottom = 'a*x_bottom^2+b*y_bottom^2+d_bottom/sqrt(1+4*a^2*x_bottom^2+4*b^2*y_bottom^2)';
        bottom_expression = {'x_bottom','y_bottom',z_bottom};
        z_top = ['a*x_top^2+b*y_top^2+d_bottom/sqrt(1+4*a^2*x_top^2+4*b^2*y_top^2)+',d_top,'/sqrt(1+4*a^2*x_top^2*(-1+2*a*d_bottom/sqrt(1+4*a^2*x_top^2+4*b^2*y_top^2)^3)^2+4*b^2*y_top^2*(-1+2*b*d_bottom/sqrt(1+4*a^2*x_top^2+4*b^2*y_top^2)^3)^2)'];
        top_expression = {'x_top','y_top',z_top};
           
        % Create the bottom surface of the Passive Corrector
        model.component('comp1').geom('geom1').create('ps1', 'ParametricSurface');
        model.component('comp1').geom('geom1').feature('ps1').label('Shim_bottom');
        model.component('comp1').geom('geom1').feature('ps1').set('selresult', true);
        model.component('comp1').geom('geom1').feature('ps1').set('parname1', 'x_bottom');
        model.component('comp1').geom('geom1').feature('ps1').set('parmin1', x_pc_bottom_min); 
        model.component('comp1').geom('geom1').feature('ps1').set('parmax1', x_pc_bottom_max);
        model.component('comp1').geom('geom1').feature('ps1').set('parname2', 'y_bottom');
        model.component('comp1').geom('geom1').feature('ps1').set('parmin2', y_pc_bottom_min);
        model.component('comp1').geom('geom1').feature('ps1').set('parmax2', y_pc_bottom_max);
        model.component('comp1').geom('geom1').feature('ps1').set('coord', bottom_expression);
        model.component('comp1').geom('geom1').feature('ps1').set('rtol', 5.0E-10);
        model.component('comp1').geom('geom1').feature('ps1').set('maxknots', 20000);
        
        % Create the top surface of the Passive Corrector
        model.component('comp1').geom('geom1').create('ps2', 'ParametricSurface');
        model.component('comp1').geom('geom1').feature('ps2').label('Shim_top');
        model.component('comp1').geom('geom1').feature('ps2').set('selresult', true);
        model.component('comp1').geom('geom1').feature('ps2').set('parname1', 'x_top');
        model.component('comp1').geom('geom1').feature('ps2').set('parmin1', x_pc_top_min);
        model.component('comp1').geom('geom1').feature('ps2').set('parmax1', x_pc_top_max);
        model.component('comp1').geom('geom1').feature('ps2').set('parname2', 'y_top');
        model.component('comp1').geom('geom1').feature('ps2').set('parmin2', y_pc_top_min);
        model.component('comp1').geom('geom1').feature('ps2').set('parmax2', y_pc_top_max);
        model.component('comp1').geom('geom1').feature('ps2').set('coord', top_expression);
        model.component('comp1').geom('geom1').feature('ps2').set('rtol', 9.0E-4);
        model.component('comp1').geom('geom1').feature('ps2').set('maxknots', 20000);
        
        % Create a solid by lofting the two surfaces
        model.component('comp1').geom('geom1').create('loft1', 'Loft');
        model.component('comp1').geom('geom1').feature('loft1').set('selresult', true);
        model.component('comp1').geom('geom1').feature('loft1').selection('profile').set({});
        model.component('comp1').geom('geom1').feature('loft1').selection('startprofile').set('ps1(1)', 1);
        model.component('comp1').geom('geom1').feature('loft1').selection('endprofile').set('ps2(1)', 1);
        model.component('comp1').geom('geom1').feature('loft1').selection('guide').set({});
        model.component('comp1').geom('geom1').run('loft1');
        
        % Add Fillet nodes
        model.component('comp1').geom('geom1').create(['fil',num2str(i+2)], 'Fillet3D');                                                
        model.component('comp1').geom('geom1').feature(['fil',num2str(i+2)]).selection('edge').set('loft1', [1 2 3 5 6 7 8 9 10 12]);    
        model.component('comp1').geom('geom1').feature(['fil',num2str(i+2)]).set('radius', '0.0005');        
        model.component('comp1').geom('geom1').feature(['fil',num2str(i+2)]).set('filletsharp', true);
        model.component('comp1').geom('geom1').run(['fil',num2str(i+2)]);
end

%% Adding materials, meshes, and physics

% Add materials
% Initialize entity array for mesh and materials selections later
entity             = linspace(2,n_mag*n_mag+1,n_mag*n_mag);
% Select all magnets (entity) for the material node. Magnet material: BMN-52
model.material('mat2').selection.set(entity);
% Select magnetic Passive Corrector material: Soft Iron (Without Losses)
if use_pc
        model.material('mat3').selection.named('geom1_loft1_dom');
end

% Add physics
% Select each Cumulative Selection nodes for the corresponding Magnetic Flux Conservation nodes  (We only need 4 magnetic flux conservation nodes if using cumulative selection, this saves memory)     
model.component('comp1').physics('mfnc').feature('mfc2').selection.named('geom1_csel1_dom');
model.component('comp1').physics('mfnc').feature('mfc3').selection.named('geom1_csel2_dom');
model.component('comp1').physics('mfnc').feature('mfc4').selection.named('geom1_csel3_dom');
model.component('comp1').physics('mfnc').feature('mfc5').selection.named('geom1_csel4_dom');
if use_pc
        model.component('comp1').physics('mfnc').feature('mfc6').selection.named('geom1_loft1_dom');
end

% Add force calculation nodes to each magnet
if calculate_force
    f = waitbar(0, 'Starting');
    for i = 1:n_mag*n_mag                               % For each magnet
        model.component('comp1').physics('mfnc').create(['fcal',num2str(i+1)], 'ForceCalculation', 3);                                       % Create a Force calculation node in physics
        model.component('comp1').physics('mfnc').feature(['fcal',num2str(i+1)]).set('TorqueRotationPoint', r_mag(i,:));                      % Set the torque rotation point as the center of the magnet 
        model.component('comp1').physics('mfnc').feature(['fcal',num2str(i+1)]).set('ForceName', ['magforce',num2str(i+1)]);                 % Give the Force of magnet a name
        model.component('comp1').physics('mfnc').feature(['fcal',num2str(i+1)]).selection.named(['geom1_blk',num2str(i+1),'_dom']);          % Select the magnet
       
        % Update waitbar
        waitbar(i/(n_mag*n_mag), f, sprintf('Adding Force Calculation nodes: %d %%', floor(i/(n_mag*n_mag)*100)));
    end
end

% Add meshes
% Select Mesh nodes for all magnets (entity), and set the mesh Fineness (1 - extremely fine, 2 - extra fine, 3 - finer, ...)
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 2);                               
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hauto', 2);      % Magnet mesh
model.component('comp1').mesh('mesh1').feature('ftet1').selection.set(entity);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hauto', 2);      % Air mesh
model.component('comp1').mesh('mesh1').feature('ftet3').feature('size1').set('hauto', 1);      % Passive Corrector mesh
if use_pc
        model.component('comp1').mesh('mesh1').feature('ftet3').selection.named('geom1_loft1_dom');
end

% Set the cutpoint in result data set 
model.result.dataset('cpt1').set('pointx', cutpoint(1));
model.result.dataset('cpt1').set('pointy', cutpoint(2));
model.result.dataset('cpt1').set('pointz', cutpoint(3));

if calculate_force
        % In result, create global evaluation nodes 
        model.result.numerical.create('gev1', 'EvalGlobal');                                   % Global evaluation for force 
        model.result.numerical.create('gev2', 'EvalGlobal');                                   % Global evaluation for torque 
        
        % In result, create 3D plot groups
        model.result.create('pg12', 'PlotGroup3D');                                            % Create 3D plot group for force
        model.result('pg12').label('Magnet Force');                             
        model.result.create('pg13', 'PlotGroup3D');                                            % Create 3D plot group for torque
        model.result('pg13').label('Magnet Torque');
        
        % Initialize the force_name, torque_name array that will be used in global evaluation and plottings
        Force_name = {};                    
        Torque_name = {};
        
        f = waitbar(0, 'Starting');
        for i = 1:n_mag*n_mag                                                                  % For each magnet
            % Concatenate the force/torque name array
           [Force_name{end+1:end+3}] =  deal(['mfnc.Forcex_magforce',num2str(i+1)], ['mfnc.Forcey_magforce',num2str(i+1)], ['mfnc.Forcez_magforce',num2str(i+1)]);
           [Torque_name{end+1:end+3}] =  deal(['mfnc.Tx_magforce',num2str(i+1)], ['mfnc.Ty_magforce',num2str(i+1)], ['mfnc.Tz_magforce',num2str(i+1)]);
            
           % Add Arrow Volume nodes in 3D plot groups
            model.result('pg12').create(['arwv',num2str(i+1)], 'ArrowVolume');                                             % Create the arrow volume node in force plot
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('expr', Force_name(end-2:end));                        % Set the expression for the force plot
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('descr', 'Electromagnetic force');                     % Set the description as Electromagnetic force
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('arrowxmethod', 'coord');                              % Set the arrow positioning  method as coordinates
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('arrowymethod', 'coord');
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('arrowzmethod', 'coord');
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('xcoord', r_mag(i,1));                                 % Set the position of the arrow (tail)
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('ycoord', r_mag(i,2));  
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('zcoord', r_mag(i,3));
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('arrowlength', 'normalized');                          % Set the arrow length as normalized
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('scaleactive', true);                                  % Enable the scale factor for arrow length
            model.result('pg12').feature(['arwv',num2str(i+1)]).set('scale', 0.0001);                                      % Set the scale factor as 0.0005
        
            model.result('pg13').create(['arwv',num2str(i+1)], 'ArrowVolume');                                             % Create the arrow volume node in torque plot                              
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('expr', Torque_name(end-2:end));                       % Set the expression for the torque plot
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('descr', 'Torque');                                    % Set the description as Torque
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('arrowxmethod', 'coord');                              % Set the arrow positioning  method as coordinates
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('arrowymethod', 'coord');
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('arrowzmethod', 'coord');
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('xcoord', r_mag(i,1));                                 % Set the position of the arrow (tail)
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('ycoord', r_mag(i,2));
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('zcoord', r_mag(i,3));
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('arrowlength', 'normalized');                          % Set the arrow length as normalized
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('scaleactive', true);                                  % Enable the scale factor for arrow length
            model.result('pg13').feature(['arwv',num2str(i+1)]).set('scale', 0.05);                                        % Set the scale factor as 0.0005
           
            % Update Waitbar
            waitbar(i/(n_mag*n_mag), f, sprintf('Evaluate and Plot the magnet force: %d %%', floor(i/(n_mag*n_mag)*100)));
        end

        % Select all force/torque components in Global Eval
        model.result.numerical('gev1').set('expr', Force_name);
        model.result.numerical('gev2').set('expr', Torque_name);
end

%% Run the study
model.study('std1').run()   

%% Postprocessing

% Evaluate the Potential value PI  at the cutpoint
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Point Evaluation 1');
model.result.numerical('pev1').set('table', 'tbl1');
model.result.numerical('pev1').setResult;
PI_level = model.result.table('tbl1').getReal;

% In each 2D plot group node, set the contour level as the evaluated PI value
model.result('pg4').feature('con1').set('levels', PI_level);
model.result('pg4').feature('con2').set('levels', PI_level);
model.result('pg5').feature('con1').set('levels', PI_level);          
model.result('pg6').feature('con1').set('levels', PI_level);
model.result('pg7').feature('con1').set('levels', PI_level);
model.result('pg8').feature('con1').set('levels', PI_level);
model.result('pg9').feature('con1').set('levels', PI_level);
model.result('pg10').feature('con1').set('levels', PI_level);
model.result('pg11').feature('iso1').set('levels', PI_level);


% Post-processing electromagnetism force and torque about the geometric center of each magnet
if calculate_force
        % Evaluate the Forces
        model.result.table.create('tbl2', 'Table');
        model.result.table('tbl2').comments('Global Evaluation 1');
        model.result.numerical('gev1').set('table', 'tbl2');
        model.result.numerical('gev1').setResult;
        model.result.table.create('tbl3', 'Table');
        model.result.table('tbl3').comments('Global Evaluation 2');
        model.result.numerical('gev2').set('table', 'tbl3');
        model.result.numerical('gev2').setResult;

        % Extract the evaluated force and torques from the model
        Force_vector = model.result.table('tbl2').getReal;
        Torque_vector = model.result.table('tbl3').getReal;
        
        % Generate X, Y coordinates using the magnet positions
        X = reshape(r_mag(:,1),[n_mag,n_mag]);
        Y = reshape(r_mag(:,2),[n_mag,n_mag]);
        
        % Form the matrix of Fx, Fy, Fz at each magnet position
        Force_vector = reshape(Force_vector,[3,n_mag*n_mag]);
        writematrix(Force_vector,"data/Force_map.txt");                             
        Fx = reshape(Force_vector(1,:),[n_mag,n_mag]);
        Fy = reshape(Force_vector(2,:),[n_mag,n_mag]);
        Fz = reshape(Force_vector(3,:),[n_mag,n_mag]);
        
        % Form the matrix of Tx, Ty, Tz at each magnet position
        Torque_vector = reshape(Torque_vector,[3,n_mag*n_mag]);
        writematrix(Force_vector,"data/Torque_map.txt");
        Tx = reshape(Torque_vector(1,:),[n_mag,n_mag]);
        Ty = reshape(Torque_vector(2,:),[n_mag,n_mag]);
        Tz = reshape(Torque_vector(3,:),[n_mag,n_mag]);
        
        % Plotting forces maps and torque maps
        figure                                                        % Figure for Fx
        surf(X, Y, Fx); 
        xlabel('X (m)','FontSize',26);                                
        ylabel('Y (m)','FontSize',26);                               
        zlabel('Fx (N)','FontSize',26);                               
        cb = colorbar;                                                % Add a colorbar
        Fx_max = max(Fx,[],"all");                                    % Find the maximum value of Fx
        Fx_min = min(Fx,[],"all");                                    % Find the minimum value of Fx 
        Fx_range = max(abs([Fx_max,Fx_min]));
        set(gca,'CLim',[-Fx_range Fx_range]);                         % Set the color scale
        ylabel(cb,'F_x (N)','FontSize',26,'Rotation',270)
        xlim([min(X,[],"all") max(X,[],"all")])
        ylim([min(Y,[],"all") max(Y,[],"all")])
        title('Electromagnet force experienced by magnets (F_x)','FontSize',26)
        axis equal
        view(2)                                                       % Set camera view from top
        
        
        figure                                                        % Figure for Fy
        surf(X, Y, Fy); 
        xlabel('X (m)','FontSize',26);                                
        ylabel('Y (m)','FontSize',26);                                
        zlabel('Fy (N)','FontSize',26);                               
        cb = colorbar;                                                % Add a colorbar 
        Fy_max =max(Fy,[],"all");                                     % Find the maximum value of Fy
        Fy_min = min(Fy,[],"all");                                    % Find the minimum value of Fy
        Fy_range = max(abs([Fy_max,Fy_min]));
        set(gca,'CLim',[-Fy_range Fy_range]);                         % Set the color scale
        ylabel(cb,'F_y (N)','FontSize',26,'Rotation',270)
        xlim([min(X,[],"all") max(X,[],"all")])
        ylim([min(Y,[],"all") max(Y,[],"all")])
        title('Electromagnet force experienced by magnets (F_y)','FontSize',26)
        axis equal
        view(2)                                                       % Set camera view from top
        
        
        figure                                                        % Figure for Fz
        surf(X, Y, Fz); 
        xlabel('X (m)','FontSize',26); 
        ylabel('Y (m)','FontSize',26); 
        zlabel('Fz (N)','FontSize',26); 
        cb = colorbar;                                                % Add a colorbar 
        Fz_max = max(Fz,[],"all");                                    % Find the maximum value of Fz
        Fz_min = min(Fz,[],"all");                                    % Find the minimum value of Fz
        Fz_range = max(abs([Fz_max,Fz_min]));
        set(gca,'CLim',[-Fz_range Fz_range]);                         % Set the color scale
        ylabel(cb,'F_z (N)','FontSize',26,'Rotation',270)
        xlim([min(X,[],"all") max(X,[],"all")])
        ylim([min(Y,[],"all") max(Y,[],"all")])
        title('Electromagnet force experienced by magnets (F_z)','FontSize',26)
        axis equal
        view(2)                                                       % Set camera view from top
        
        
        figure;                                                       % Figure for Tx
        surf(X, Y, Tx); 
        xlabel('X (m)','FontSize',26); 
        ylabel('Y (m)','FontSize',26); 
        zlabel('Tx (N*m)','FontSize',26); 
        cb = colorbar;                                                % Add a colorbar 
        Tx_max = max(Tx,[],"all");                                    % Find the maximum value of Tx
        Tx_min = min(Tx,[],"all");                                    % Find the maximum value of Tx
        Tx_range = max(abs([Tx_max,Tx_min]));
        set(gca,'CLim',[-Tx_range Tx_range]);                         % Set the color scale
        ylabel(cb,'T_x (N*m)','FontSize',26,'Rotation',270)
        xlim([min(X,[],"all") max(X,[],"all")])
        ylim([min(Y,[],"all") max(Y,[],"all")])
        title('Torque about magnet center (T_x)','FontSize',26)
        axis equal
        view(2)                                                       % Set camera view from top
        
        
        figure;                                                       % Figure for Ty
        surf(X, Y, Ty); 
        xlabel('X (m)','FontSize',26);     
        ylabel('Y (m)','FontSize',26);      
        zlabel('Ty (N*m)','FontSize',26); 
        cb = colorbar;                                                % Add a colorbar
        Ty_max = max(Ty,[],"all");                                    % Find the maximum value of Ty
        Ty_min = min(Ty,[],"all");                                    % Find the maximum value of Ty
        Ty_range = max(abs([Ty_max,Ty_min]));
        set(gca,'CLim',[-Ty_range Ty_range]);                         % Set the color scale
        ylabel(cb,'T_y (N*m)','FontSize',26,'Rotation',270)
        xlim([min(X,[],"all") max(X,[],"all")])
        ylim([min(Y,[],"all") max(Y,[],"all")])
        title('Torque about magnet center (T_y)','FontSize',26)
        axis equal
        view(2)                                                       % Set camera view from top
        
        figure;                                                       % Figure for Tz
        surf(X, Y, Tz); 
        xlabel('X (m)','FontSize',26);      
        ylabel('Y (m)','FontSize',26);      
        zlabel('Tz (N*m)','FontSize',26);
        cb = colorbar;                                                % Add a colorbar
        Tz_max = max(Tz,[],"all");                                    % Find the maximum value of Tz
        Tz_min = min(Tz,[],"all");                                    % Find the maximum value of Tz
        Tz_range = max(abs([Tz_max,Tz_min]));
        set(gca,'CLim',[-Tz_range Tz_range]);                         % Set the color scale
        ylabel(cb,'T_z (N*m)','FontSize',26,'Rotation',270)
        xlim([min(X,[],"all") max(X,[],"all")])
        ylim([min(Y,[],"all") max(Y,[],"all")])
        title('Torque about magnet center (T_z)','FontSize',26)
        axis equal
        view(2)                                                       % Set camera view from top
    
        % Plot Force map
        figure
        pd12 = mphplot(model,"pg12");
        title('Electromagnet force experienced by magnets')
        
        % Plot Torque map
        figure
        pd13 = mphplot(model,"pg13");
        title('Torque about magnet center')
end


% Post-processing Equipotential lines at a centain distance above the Halbach Array
Equipotential_surface = 1;                                                    %  Set to '1' to plot the contour plots, '0' to not plot
if Equipotential_surface
        figure
        pd4 = mphplot(model,"pg4");
        figure
        pd5 = mphplot(model,"pg5");
        figure
        pd6 = mphplot(model,"pg6");
        figure
        pd7 = mphplot(model,"pg7");
        figure
        pd8 = mphplot(model,"pg8");
        figure
        pd9 = mphplot(model,"pg9");
        figure
        pd10 = mphplot(model,"pg10");
        
        % Define the Ideal surface
        x_ideal = linspace(-0.25,0.25);                                     % [m]
        r_sphere = 2;                                                       % [m]
        focalLength = 1;                                                    % [m]
        z_ideal = -sqrt(r_sphere^2-x_ideal.^2)+r_sphere+h_ff2hb+delta_ff;   % [m]
        
        % Define a invisible scan sphere radius centered at the focal point
        ScanSphere_radius = 1.98;                                           % Scan Sphere surface [m]
        
        % Extract the data points from the plot group contour using a Scan
        % Sphere for each cut plane (0 deg, 30 deg, 60 deg, 90 deg, 120 deg, 150 deg, 180 deg)
        p_0deg = pd5{1,3}{1,1}.p;
        P_0deg = extractpoints(p_0deg,[0;2*focalLength],ScanSphere_radius);
        
        p_30deg = pd6{1,3}{1,1}.p;
        P_30deg = extractpoints(p_30deg,[0;2*focalLength],ScanSphere_radius);
        
        p_60deg = pd7{1,3}{1,1}.p;
        P_60deg = extractpoints(p_60deg,[0;2*focalLength],ScanSphere_radius);
        
        p_90deg = pd8{1,3}{1,1}.p;
        P_90deg = extractpoints(p_90deg,[0;2*focalLength],ScanSphere_radius);
        
        p_120deg = pd9{1,3}{1,1}.p;
        P_120deg = extractpoints(p_120deg,[0;2*focalLength],ScanSphere_radius);
        
        p_150deg = pd10{1,3}{1,1}.p;
        P_150deg = extractpoints(p_150deg,[0;2*focalLength],ScanSphere_radius);
        
        h_compensation = interp1(P_0deg(1,size(P_0deg,2)/3:size(P_0deg,2)*2/3),P_0deg(2,size(P_0deg,2)/3:size(P_0deg,2)*2/3),0,"linear")-h_ff2hb-delta_ff;
        
        % Plot the equipotential surfaces together with the ideal surface
        figure
        plot(x_ideal,z_ideal,'r--',P_0deg(1,:),P_0deg(2,:)-h_compensation,P_30deg(1,:),P_30deg(2,:)-h_compensation,P_60deg(1,:),P_60deg(2,:)-h_compensation,P_90deg(1,:),P_90deg(2,:)-h_compensation,P_120deg(1,:),P_120deg(2,:)-h_compensation,P_150deg(1,:),P_150deg(2,:)-h_compensation)
        set(0, 'DefaultLineLineWidth', 1.5);
        grid on
        xlabel('X (m)','FontSize',20)
        ylabel('Y (m)','FontSize',20)
        title('Equipotential lines at z = 11 mm, in symmetric cut-planes','FontSize',24)
        legend('Ideal surface','0 deg (yz-plane)','30 deg','60 deg','90 deg (xz-plane)','120 deg','150 deg',fontsize = 16)
        ylim([0 0.03])

        figure
        plot(x_ideal,z_ideal,'r--',P_0deg(1,:),P_0deg(2,:)-h_compensation,P_90deg(1,:),P_90deg(2,:)-h_compensation)
        set(0, 'DefaultLineLineWidth', 1.5);
        grid on
        xlabel('X (m)','FontSize',20)
        ylabel('Y (m)','FontSize',20)
        title('Equipotential lines at z = 11 mm, in symmetric cut-planes','FontSize',24)
        legend('Ideal surface','0 deg (yz-plane)','90 deg (xz-plane)',fontsize = 16)
        ylim([0 0.03])
end

% General plot setting
set(findobj('Type','axes'),'FontName','Arial','FontSize',28,'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter','latex','TickDir','out');
set(findobj('Type','Legend'),'Interpreter','latex','box','off');

%% Save and open the model
% Save the model
mphsave(model,'HalbachArraySimulation.mph')

% Open the model in COMSOL GUI
mphlaunch(model)

toc