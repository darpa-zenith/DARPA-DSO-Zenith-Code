function [bx,by,bz] = magloop(x,y,z,a,h)

% a: coil radius
% h: coil vertical position?

r = sqrt(x.^2+y.^2)+eps; % distance to (x,y)

u0 = 4*pi*1e-7;
k = sqrt(4*a*r./((r+a).^2+(z-h).^2));
[K,E] = ellipke(k.^2);

b0 = u0*k./(4*pi*sqrt(a*r.^3));

br = -b0.*(z-h).*(K-E.*(2-k.^2)./(2*(1-k.^2)));

bx = br.*x./r;
by = br.*y./r;

bz = b0.*r.*(K+E.*(k.^2.*(r+a)-2*r)./(2*r.*(1-k.^2)));

