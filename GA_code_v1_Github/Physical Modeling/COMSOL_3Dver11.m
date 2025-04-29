clear all
clc

A=xlsread('0.5m Full Parabola 4mm equal slew 1 deg per sec END.csv');
A=[A(:,1) A(:,2)]; 
x0=A(:,1);
z0=A(:,2);

% The theta defines 12 different rotational angles to revolve
% the surface lineout of Mercury from the file in line 4.
theta=linspace(0.2,2*pi,10*pi);%degrees

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

% 7
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
% 13
x13=x0*cos(theta(13));
y13=x0*sin(theta(13));
 
xtest=cat(1,xtest,x13);
ytest=cat(1,ytest,y13);
ztest=cat(1,ztest,z0);

%% 14
x14=x0*cos(theta(14));
y14=x0*sin(theta(14));
 
xtest=cat(1,xtest,x14);
ytest=cat(1,ytest,y14);
ztest=cat(1,ztest,z0);

%% 15
x15=x0*cos(theta(15));
y15=x0*sin(theta(15));
 
xtest=cat(1,xtest,x15);
ytest=cat(1,ytest,y15);
ztest=cat(1,ztest,z0);

%% 16
x16=x0*cos(theta(16));
y16=x0*sin(theta(16));
 
xtest=cat(1,xtest,x16);
ytest=cat(1,ytest,y16);
ztest=cat(1,ztest,z0);

%% 17
x17=x0*cos(theta(17));
y17=x0*sin(theta(17));
 
xtest=cat(1,xtest,x17);
ytest=cat(1,ytest,y17);
ztest=cat(1,ztest,z0);

%% 18
x18=x0*cos(theta(18));
y18=x0*sin(theta(18));
 
xtest=cat(1,xtest,x18);
ytest=cat(1,ytest,y18);
ztest=cat(1,ztest,z0);

%% 19
x19=x0*cos(theta(19));
y19=x0*sin(theta(19));

xtest=cat(1,xtest,x19);
ytest=cat(1,ytest,y19);
ztest=cat(1,ztest,z0);
%% 20

x20=x0*cos(theta(20));
y20=x0*sin(theta(20));
 
xtest=cat(1,xtest,x20);
ytest=cat(1,ytest,y20);
ztest=cat(1,ztest,z0);

%% 21
x21=x0*cos(theta(21));
y21=x0*sin(theta(21));
 
xtest=cat(1,xtest,x21);
ytest=cat(1,ytest,y21);
ztest=cat(1,ztest,z0);
%% 22
x22=x0*cos(theta(22));
y22=x0*sin(theta(22));
 
xtest=cat(1,xtest,x22);
ytest=cat(1,ytest,y22);
ztest=cat(1,ztest,z0);

%% 23
x23=x0*cos(theta(23));
y23=x0*sin(theta(23));
 
xtest=cat(1,xtest,x23);
ytest=cat(1,ytest,y23);
ztest=cat(1,ztest,z0);

%% 24
x24=x0*cos(theta(24));
y24=x0*sin(theta(24));
 
xtest=cat(1,xtest,x24);
ytest=cat(1,ytest,y24);
ztest=cat(1,ztest,z0);

%% 25
x25=x0*cos(theta(25));
y25=x0*sin(theta(25));
 
xtest=cat(1,xtest,x25);
ytest=cat(1,ytest,y25);
ztest=cat(1,ztest,z0);

%% 26
x26=x0*cos(theta(26));
y26=x0*sin(theta(26));
 
xtest=cat(1,xtest,x26);
ytest=cat(1,ytest,y26);
ztest=cat(1,ztest,z0);

%% 27
x27=x0*cos(theta(27));
y27=x0*sin(theta(27));
 
xtest=cat(1,xtest,x27);
ytest=cat(1,ytest,y27);
ztest=cat(1,ztest,z0);

%% 28
x28=x0*cos(theta(28));
y28=x0*sin(theta(28));
 
xtest=cat(1,xtest,x28);
ytest=cat(1,ytest,y28);
ztest=cat(1,ztest,z0);

%% 29
x29=x0*cos(theta(29));
y29=x0*sin(theta(29));
 
xtest=cat(1,xtest,x29);
ytest=cat(1,ytest,y29);
ztest=cat(1,ztest,z0);

%%----------------------------------------
%% 30
x30=x0*cos(theta(30));
y30=x0*sin(theta(30));
 
xtest=cat(1,xtest,x30);
ytest=cat(1,ytest,y30);
ztest=cat(1,ztest,z0);

% 31
x31=x0*cos(theta(31));
y31=x0*sin(theta(31));
 
xtest=cat(1,xtest,x31);
ytest=cat(1,ytest,y31);
ztest=cat(1,ztest,z0);

%% 32
% x32=x0*cos(theta(32));
% y32=x0*sin(theta(32));
%  
% xtest=cat(1,xtest,x32);
% ytest=cat(1,ytest,y32);
% ztest=cat(1,ztest,z0);
% 
% %% 33
% x33=x0*cos(theta(33));
% y33=x0*sin(theta(33));
%  
% xtest=cat(1,xtest,x33);
% ytest=cat(1,ytest,y33);
% ztest=cat(1,ztest,z0);
% 
% %% 34
% x34=x0*cos(theta(34));
% y34=x0*sin(theta(34));
%  
% xtest=cat(1,xtest,x34);
% ytest=cat(1,ytest,y34);
% ztest=cat(1,ztest,z0);
% 
% %% 35
% x35=x0*cos(theta(35));
% y35=x0*sin(theta(35));
%  
% xtest=cat(1,xtest,x35);
% ytest=cat(1,ytest,y35);
% ztest=cat(1,ztest,z0);
%%----------------------------------------

% xlin = linspace(min(xtest), max(xtest), 100);
% ylin = linspace(min(ytest), max(ytest), 100);
% [X,Y] = meshgrid(xlin, ylin);
% 
% Z = griddata(xtest,ytest,ztest,X,Y,'v4');
% mesh(X,Y,Z)
% axis tight; hold on
% plot3(xtest,ytest,ztest,'.','MarkerSize',15) 
% title('Residual Surface Error')
% xlabel('mm');
% ylabel('mm');
% zlabel('mm');
%%----------------------------------------
A=[xtest ytest ztest]; 
r=sqrt(xtest.^2+ytest.^2);

xtest=xtest;
ytest=ytest;
ztest=ztest;

format long
zz=0.7503 + 0.0003334*xtest.^2 + 0.0003334*ytest.^2;

q=ztest-zz;
A=[xtest ytest q];

B=[ round(A(:,1)) round(A(:,2))  A(:,3)];
x=B(:,1);
y=B(:,2);
z=B(:,3);
figure;
 
%%----------------------------------------

xlin = linspace(min(x), max(x), 600);
ylin = linspace(min(y), max(y), 600);
[X,Y] = meshgrid(xlin, ylin);

Z = griddata(x,y,z,X,Y,'v4');
mesh(X,Y,Z)
axis tight; hold on

%% NEW
[m,n]=size(X);
[X,Y]=meshgrid(1:n,1:m);
P=[X(:) Y(:) Z(:)*1000]

%% NEW
zeta=accumarray(P(:,[2 1]),P(:,3));
writematrix(zeta);
type 'zeta.txt';

%Data = readtable('ExcelFile.csv');
%writetable(Data, 'textfile.txt');


% Plot of data and control of range along z-axis
title('Residual Surface Error')
xlabel('mm');
ylabel('mm');
zlabel('mm');
zlim([-0.1 0.1])
