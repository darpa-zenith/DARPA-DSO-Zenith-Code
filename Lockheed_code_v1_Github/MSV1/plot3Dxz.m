function plot3Dxz(Hval, xval, yval, zval)
%PLOT3DXZ Plot values of Hval
%   Plot in X,Z with angle fixed at zero

size(Hval)
size(xval)
size(yval)
size(zval)

numx = size(Hval, 1);
numy = size(Hval, 2);
numz = size(Hval, 3);

slen = 10;

x2 = zeros(slen, 1);
y2 = zeros(slen, 1);
z2 = zeros(slen, 1);
    
xx = rand(slen, 2);
x2 = 0.25 * xx(:, 1);
y2 = zeros(slen, 1);
z2 = 0.01 * xx(:, 2) + 0.0075;

figure
hold on

for k = 1:numz
    xi = xval(:,1,k);
    zi = zval(:,1,k);
    Hi = Hval(:,1,k);
    plot3(xi, zi, Hi, 'k-');
end

if 1
for j = 1:10:numx
    for k = 1:numz
        xj(k)= xval(j,1,k);
        zj(k) = zval(j,1,k);
        Hj(k) = Hval(j,1,k);
    end
    plot3(xj, zj, Hj, 'r-');
end
end


pts2 = interpXYZ(xval, yval, zval, Hval, x2, y2, z2);

plot3(x2(:,1), z2(:,1), pts2, 'go');
plot3(x2(:,1), z2(:,1), pts2, 'kx');

xlabel('radial position (m)')
ylabel('height (m)')
zlabel('H\-mag')

% pts
% pts2
% pts2 - pts

end