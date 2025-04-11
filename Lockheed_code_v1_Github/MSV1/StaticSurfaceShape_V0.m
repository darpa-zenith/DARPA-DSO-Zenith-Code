% Static Surface Steps:
% - Set up parabolic and cartesian grids   
% - Compute body force loading on grids
% - Read in magnetic potential, interpolate onto parabolic grid
% - Compute laplace matrix coefficients
% - Initialize surface height
% - loop
%    - compute source terms on surface
%    - solve for new surface
%    - compute change in surface magnitude
%    - if change is "small", terminate loop
% - write out surface shape
% - plot surface shape
