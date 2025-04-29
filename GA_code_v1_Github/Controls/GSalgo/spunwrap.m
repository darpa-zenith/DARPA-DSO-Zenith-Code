function [out, lsout, bpout] = spunwrap(in, recon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the phase unwrapping for an input wrapped phase
% map. It computes a least squares solution to the first phase difference
% of the phase points using a set of predefined matrices. These matrices
% can be readily created using the mksprecon.m function. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [out, lsout, bpout] = spunwrap(in, recon)
% -- Inputs --
%   in = [M, N] Double : The array of the wrapped phase that is to be
%       unwrapped.
%   recon = Structure : Contains the matrices and index mapping used in
%       the least squares reconstruction of the phase. See the 
%       mksprecon function to see how these fields are computed. This 
%       structure contains the following fields:
%           g : Matrix representation of X and Y differences.
%           gt : g transpose.
%           p : symamd(g'*g), for a symmetric positive definite matrix 
%               S = G'*G, returns the permutation vector p such 
%               that S(p,p) tends to have a sparser Cholesky factor
%               than S.
%           r : chol(S), where S = gtg(p,p) (gtg = G'*G)
%               r is the sparse upper triangular matrix s.t. S =
%               r'*r.
%           rt : r transpose.
%           maskx : Mask for x gradients.
%           masky : Mask for y gradients.
%           maska : Mask for all points where the phase is unwrapped.
%           indx  : Indices for x gradients.
%           indy  : Indices for y gradients.
%           inda  : Indices for all unwrapped points.
% 
% -- Outputs --
%   out = [M, N] Double : The total unwrapped phase map. This is the
%       sum of the least squares and the branch point phase.
%   lsout = [M, N] Double : The least squares (ls) portion of the 
%       unwrapped phase.
%   bpout = [M, N] Double : The branch point (bp) portion of the
%       unwrapped phase.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Condition inputs
nx = size(in);
if nargin < 2
    recon = mksprecon(find(in),nx);
end

% Compute Gradients
sx = angle(exp(1i*in(:,2:end)).*exp(-1i*in(:,1:end-1)));
sy = angle(exp(1i*in(2:end,:)).*exp(-1i*in(1:end-1,:)));
st = [sx(recon.indx); sy(recon.indy)];
% Initialize unwrapped phase arrays
out = zeros(nx(1),nx(2));
lsout = zeros(nx(1),nx(2));
bpout = zeros(nx(1),nx(2));
% Compute least square phase
lsout(recon.inda) = sprecon(st,recon.gt,recon.r,recon.rt,recon.p);
% Compute branch point phase ("hidden phase")
bpout(recon.inda) = angle(exp(1i*in(recon.inda)).*exp(-1i*lsout(recon.inda)));
% Compute total phase
out(recon.inda) = lsout(recon.inda) + bpout(recon.inda);

end
