function a = sprecon(s, gt, r, rt, p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the least squares solution based on a Cholesky
% decomposition of the gradient solution matrix and a vector of X and Y
% phase gradient vectors. These matrices can be readily computed using the
% mksprecon for a particular aperture geometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a = sprecon(s, gt, r, rt, p)
% -- Inputs --
%   s = [M, 1] Double : The vector of X and Y phase differences.
%   gt = Array : Transpose of the gradient solution matrix g.
%   r = The Cholesky decomposition of the least squares solution to the
%       gradient solution g'*g.
%   rt = The transpose of the matrix r.
%   p = Set of indexed points to allow the Cholesky decomposition to be
%       more sparse.
% 
% -- Outputs --
%   a = Double : The solution of the least squares gradient matrix
%       g*a = s. Rather, this is the vector of the unwrapped phase
%       reconstruction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute g'*s to solve (g'*g)*a = g'*s
gts=gt*s;
% Index in same way as Cholesky decomposition for sparser representation
y = gts(p);
% Solve for a
yp=rt\y;
ap=r\yp;
% Revert index ordering to match initial spatial indexing
a(p)=ap;
a=a';

end
