function plot3Dxy(Hval, xval, yval, zval)
%PLOT3DXY Plot values of Hval
%   generates slen random locations in x and z
%   plots the surface, interpolates the values
%   at the random locations and plots them on 
%   the surface.  
%   The z value is fixed and can be modified by 
%   adjusting z2.

%   Assume a 2-D surface, in x and z 

size(Hval)
size(xval)
size(yval)
size(zval)

numx = size(Hval, 1);
numy = size(Hval, 2);
numz = size(Hval, 3);

slen = 10;
numSamp = 0;

x2 = zeros(slen, 1);
y2 = zeros(slen, 1);
z2 = zeros(slen, 1);

while numSamp < slen
    
    xx = rand(1, 2, 3);
    x1 = 0.5 * xx(1, 1) - 0.25;
    y1 = 0.5 * xx(1, 2) - 0.25;
    if norm([x1, y1])  <= 0.25
        numSamp = numSamp + 1;
        x2(numSamp) = x1;
        y2(numSamp) = y1;
        % z2(numSamp) = 0.01 * xx(1, 3) + 0.0075;
        z2(numSamp) = 0.0174;
    end

end

figure
hold on

for k = 1:numy
    plot3(xval(:,k,1), yval(:,k,1), Hval(:,k,21), '-');
end

for j = 1:10:numx
    plot3(xval(j,:,1), yval(j,:,1), Hval(j,:,21), '-');
end

pts2 = interpXYZ(xval, yval, zval, Hval, x2, y2, z2);

plot3(x2(:,1), y2(:,1), pts2, 'cx');
plot3(x2(:,1), y2(:,1), pts2, 'go');

xlabel('x position (m)');
ylabel('y position (m)');
zlabel('H\_mag');


end