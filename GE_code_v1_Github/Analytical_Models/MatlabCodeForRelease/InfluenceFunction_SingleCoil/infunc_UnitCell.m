clc;clear all;close all;

[X,Y] = meshgrid(linspace(-0.08,0.08,100));

a = 2*[0.0175:-0.00081:0.0039]; % coil geometry - diameters of wires
h = [-0.066:0.00081:-0.032]; % coil geometry - heights of wires
current = 0.1;
offset = -0.035;
h = h + offset;
BX = 0; BY = 0; BZ = 0;
for i = 1:length(h)
    for j = 1:length(a)
       [bx,by,bz] = magloop(X,Y,0,a(j),h(i));
       BX = BX + bx;
       BY = BY + by;
       BZ = BZ + bz;
    end
end

ur = 2.7;           
rho = 1210;        
u0 = 4*pi*1e-7; 
g = 9.8;
k = (ur-1)/(u0*ur*rho*g);
w = k*(ur*(BX).^2+ur*(BY).^2+(BZ).^2);

figure;
plot(X(50,:), current*BZ(50,:)*10000) % unit Gauss
xlabel('coordinates (m)')
ylabel('Magnetic Field (Gauss)')

figure;
plot(X(50,:), w(50,:)*1e6) % unit Gauss
xlabel('coordinates (m)')
ylabel('Deformation (um)')