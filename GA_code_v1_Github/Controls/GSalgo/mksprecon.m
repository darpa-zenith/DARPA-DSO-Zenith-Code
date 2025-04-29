function recon = mksprecon(indap, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes a set of matrices used to unwrap a spatial array
% of phase information. These matrices are constructed to efficiently
% compute the least squares fit of the slopes, assuming continuous phase,
% using a set of sparse matrices. These matrices are stored in a structure
% to be used in tandem with the spunwrap function to compute the unwrapped
% phase pattern.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recon = mksprecon(indap, m)
% -- Inputs --
%   indap = [M, 1] Integer : The index locations of the array that are
%       to be unwrapped.
%   m = Integer : The total size of the 2D array over which the phase 
%       is unwrapped.
% 
% -- Outputs --
%   recon = Structure : Contains the matrices and index mapping used in
%       the least squares reconstruction of the phase. This structure
%       contains the following fields:
%           g : Matrix representation of X and Y differences.
%           gt : g transpose.
%           p : symamd(g'*g), for a symmetric positive definite matrix 
%                 S = G'*G, returns the permutation vector p such that 
%                 S(p,p) tends to have a sparser Cholesky factor 
%                 than S.
%           r : chol(S), where S = gtg(p,p) (gtg = G'*G)
%                 r is the sparse upper triangular matrix s.t. S =
%                 r'*r.
%           rt : r transpose.
%           maskx : Mask for x gradients.
%           masky : Mask for y gradients.
%           maska : Mask for all points where the phase is unwrapped.
%           indx  : Indices for x gradients.
%           indy  : Indices for y gradients.
%           inda  : Indices for all unwrapped points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(numel(m) ==1)
    m = [m, m];
end
% Compute index ordering in geometric 2D area
ap = zeros(m);
ap(indap) = (1:length(indap));

% Determine points used for difference in X
ap1x = ap(:,1:m(2)-1);
ap2x=ap(:,2:m(2));
% Determine points used for difference in Y
ap1y=ap(1:m(1)-1,:);
ap2y=ap(2:m(1),:);

% Index points that have a first difference in defined indexed map
xm=find(ap1x > 0 & ap2x > 0);
ym=find(ap1y > 0 & ap2y > 0);

% Define total number of points with differences
lxm = length(xm);
lym = length(ym);

% Create sparse array defining matrix multiply represenation of X and Y
%   difference
j = [ap1x(xm)
     ap2x(xm)
     ap1y(ym)
     ap2y(ym)];
i = [(1:lxm)'
     (1:lxm)'
     (lxm+1:lym+lxm)'
     (lxm+1:lym+lxm)'];
s = [-ones(lxm,1)
      ones(lxm,1)
     -ones(lym,1)
      ones(lym,1)];
g = sparse(i,j,s);
gt = g';
% Compute LS matrix solution (g'*g)
gtg = gt*g;
% Determine reording indices of g'*g for sparse cholesky decomposition
p = symamd(gtg);
% Add identity matrix for stability in solution
gtgperm = gtg(p,p) + 1.0e-10*speye(size(gtg));
% Compute Cholesky decomposition of LS matrix solution g'*g
r = chol(gtgperm);
rt = r';

% Define masking function for X difference
maskx = zeros(m(1),m(2)-1);
maskx(xm) = 1;
% Define masking function for Y difference
masky = zeros(m(1)-1,m(2));
masky(ym) = 1;
% Define total masking function
maska = zeros(m(1),m(2));
maska(indap) = 1;

% Store matrices in output structure
recon.g = g;
recon.gt = gt;
recon.p = p;
recon.r = r;
recon.rt = rt;
recon.maskx = maskx;
recon.masky = masky;
recon.maska = maska;
recon.indx = find(maskx);
recon.indy = find(masky);
recon.inda = find(maska);

end
