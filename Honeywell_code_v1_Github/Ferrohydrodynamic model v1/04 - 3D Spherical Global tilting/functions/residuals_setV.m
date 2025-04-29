function Vres = residuals_setV(z0, params)
% Computes the volume of the interface whose 0,0 height is at z0 minus the
% target volume
%
% Vres = volumeinterface(z0, params)
% 
% Inputs:
%           z0      Target height at (0,0)
%           params  Parameters structure
%
% Outputs:
%           Vres    Interface volume residual

% Expand parameters structure
z_dish    = params.z_dish;
x_t       = params.x_t;
y_t       = params.y_t;
A         = params.A;
dz        = params.dz;
V0        = params.V0;
R         = params.R;

% Tilted interface starting at z = z0
opt       = optimoptions('fsolve','Display','none');
z_t       = fsolve(@(Z)residuals_setPI(Z, z0, params), -sqrt(R^2 - x_t.^2 - y_t.^2),opt)  - dz;

% Remove outer region
outregion      = sqrt(x_t.^2+y_t.^2) > A/2;
z_t(outregion) = z_dish(outregion);

% Compute Volume
Vres         = 1e7*(quad2d(@(x,y) interp2(x_t',y_t',z_t'-z_dish',x,y,'linear'),...
            min(min(x_t)),max(max(x_t)),min(min(y_t)),max(max(y_t)),'AbsTol',1e-8) - V0);
