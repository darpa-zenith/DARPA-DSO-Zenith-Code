function Plot_magnets(theta,w,h,r0,cover)
%  This function plot the magnets shape with input:
%  theta (Magnets orientations with respect to vertical axis) [rad] 
%
%  Plot_magnets(theta,w,h,r0,cover)
%
%  Input  : 
%           theta (Magnets orientations with respect to
%                  vertical axis)                       [rad] 
%           w     (Widths of magnets)                   [m]
%           h     (Heights of magnets)                  [m]
%           r0    (Position vectors of each magnets)    [m]
%`          cover (whether to plot (1) or not (0))
%
%  Output : plot representation of magnets over the existing plot


R           = [cos(theta) -sin(theta);sin(theta) cos(theta)];
direc       = R*[0;1];
r1          = R*([-w/2;0]+[0;h/2])+r0';
r2          = R*([-w/2;0]-[0;h/2])+r0';
r3          = R*([w/2;0]+[0;h/2])+r0';
r4          = R*([w/2;0]-[0;h/2])+r0';

xcord       =[r1(1),r2(1),r4(1),r3(1),r1(1)];
ycord       =[r1(2),r2(2),r4(2),r3(2),r1(2)];
plot(xcord,ycord,'LineWidth',1.5,'Color','k')
plot(r0(1),r0(2),'+m','LineWidth',1.5)

c           = [112,128,144];
c           = c./norm(c);
if cover
    patch(xcord(1:end-1),ycord(1:end-1),c)
end

quiver(r0(1),r0(2),direc(1),direc(2), 0.5*min([w,h]), 'k', 'LineWidth', 1.5,'MaxHeadSize',1.5)
end

