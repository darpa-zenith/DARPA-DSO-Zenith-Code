function [bx,by,bz] = magloop(x,y,z,a,h)
% Function that computes the magnetic field produce by a current loop at any point in space
%
% INPUTS:
%  x, y, z: x, y and z coordinates at which the magnetic field is calculated (scalar or vector)
%        a: Coil radius (scalar)
%        h: z offset position of the coil (scalar)
% OUTPUTS:
% bx, by, bz: Values of the magnetic field components in gauss units (scalar or vector)

% Equations adapted from the following sources:
% https://netdenizen.com/emagnet/offaxis/off_axis_loop.html
% "SOME USEFUL INFORMATION FOR THE DESIGN OF AIR-CORE SOLENOIDS," D.Bruce Montgomery and J. Terrell, Nov. 1961 under contract from Air Force AF19(604)-7344 
% "STATIC AND DYNAMIC ELECTRICITY," W.R. Smythe, McGraw-Hill, New York, 1950, p. 266
 
r = sqrt(x.^2+y.^2)+eps;
u0 = 4*pi*1e-7;
k = sqrt(4*a*r./((r+a).^2+(z-h).^2));
[K,E] = ellipke(k.^2);
b0 = u0*k./(4*pi*sqrt(a*r.^3));
br = -b0.*(z-h).*(K-E.*(2-k.^2)./(2*(1-k.^2)));
bx = br.*x./r;
by = br.*y./r;
bz = b0.*r.*(K+E.*(k.^2.*(r+a)-2*r)./(2*r.*(1-k.^2)));
% end
