clear all
clc

A=xlsread('Parabola 2mm Slewing end Center Tub nonzero deg.csv');
A=[A(:,1) A(:,2)]; 
x0=A(:,1);
z0=A(:,2);

% The theta defines 12 different rotational angles to revolve
% the surface lineout of Mercury from the file in line 4.
theta=linspace(0.6283,2*pi,4*pi);%degrees

% The data is then appended in two columns for x and y using 
% multiple instances of lines 22-27. Further code refinement 
% would use a for loop.
x0=x0;
y0=x0.*0;

xtest=x0;
ytest=y0;
ztest=z0;

x1=x0*cos(theta(1));
y1=x0*sin(theta(1));

xtest=cat(1,x0,x1);
ytest=cat(1,y0,y1);
ztest=cat(1,z0,z0);

%% 2
x2=x0*cos(theta(2));
y2=x0*sin(theta(2));

xtest=cat(1,xtest,x2);
ytest=cat(1,ytest,y2);
ztest=cat(1,ztest,z0);
%% 3

x3=x0*cos(theta(3));
y3=x0*sin(theta(3));
 
xtest=cat(1,xtest,x3);
ytest=cat(1,ytest,y3);
ztest=cat(1,ztest,z0);

%% 4
x4=x0*cos(theta(4));
y4=x0*sin(theta(4));
 
xtest=cat(1,xtest,x4);
ytest=cat(1,ytest,y4);
ztest=cat(1,ztest,z0);
%% 5
x5=x0*cos(theta(5));
y5=x0*sin(theta(5));
 
xtest=cat(1,xtest,x5);
ytest=cat(1,ytest,y5);
ztest=cat(1,ztest,z0);

%% 6
x6=x0*cos(theta(6));
y6=x0*sin(theta(6));
 
xtest=cat(1,xtest,x6);
ytest=cat(1,ytest,y6);
ztest=cat(1,ztest,z0);

%% 7
x7=x0*cos(theta(7));
y7=x0*sin(theta(7));
 
xtest=cat(1,xtest,x7);
ytest=cat(1,ytest,y7);
ztest=cat(1,ztest,z0);

%% 8
x8=x0*cos(theta(8));
y8=x0*sin(theta(8));
 
xtest=cat(1,xtest,x8);
ytest=cat(1,ytest,y8);
ztest=cat(1,ztest,z0);

%% 9
x9=x0*cos(theta(9));
y9=x0*sin(theta(9));
 
xtest=cat(1,xtest,x9);
ytest=cat(1,ytest,y9);
ztest=cat(1,ztest,z0);

%% 10
x10=x0*cos(theta(10));
y10=x0*sin(theta(10));
 
xtest=cat(1,xtest,x10);
ytest=cat(1,ytest,y10);
ztest=cat(1,ztest,z0);

%% 11
x11=x0*cos(theta(11));
y11=x0*sin(theta(11));
 
xtest=cat(1,xtest,x11);
ytest=cat(1,ytest,y11);
ztest=cat(1,ztest,z0);

%% 12
x12=x0*cos(theta(12));
y12=x0*sin(theta(12));
 
xtest=cat(1,xtest,x12);
ytest=cat(1,ytest,y12);
ztest=cat(1,ztest,z0);

%%----------------------------------------

xlin = linspace(min(xtest), max(xtest), 100);
ylin = linspace(min(ytest), max(ytest), 100);
[X,Y] = meshgrid(xlin, ylin);

Z = griddata(xtest,ytest,ztest,X,Y,'v4');
mesh(X,Y,Z)
axis tight; hold on
plot3(xtest,ytest,ztest,'.','MarkerSize',15) 
title('Residual Surface Error')
xlabel('mm');
ylabel('mm');
zlabel('mm');
%%----------------------------------------
A=[xtest ytest ztest]; 
r=sqrt(xtest.^2+ytest.^2);

xtest=xtest;
ytest=ytest;
ztest=ztest;

format long
zz=1.709 + 0.001988*xtest.^2 + 0.001988*ytest.^2;

q=ztest-zz;
A=[xtest ytest q];

B=[ round(A(:,1)) round(A(:,2))  A(:,3)];
x=B(:,1);
y=B(:,2);
z=B(:,3);
figure;

%%----------------------------------------

xlin = linspace(min(x), max(x), 100);
ylin = linspace(min(y), max(y), 100);
[X,Y] = meshgrid(xlin, ylin);

Z = griddata(x,y,z,X,Y,'v4');
mesh(X,Y,Z)
axis tight; hold on

% Plot of data and control of range along z-axis
title('Residual Surface Error')
xlabel('mm');
ylabel('mm');
zlabel('mm');
zlim([-0.035 0.05])
