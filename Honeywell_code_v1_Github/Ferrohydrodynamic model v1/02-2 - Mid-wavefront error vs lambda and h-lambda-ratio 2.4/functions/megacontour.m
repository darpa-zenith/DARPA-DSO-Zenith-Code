function megacontour(x, y, Z, c0, c1, color)
%MEGACONTOUR This function plots a contour inequality by plotting a solid line denoting the equality, and a dashed line denoting the side of the equality you
%want to be on.

% bpoint = megacontour(x, y, Z, c0, c1, color)
%
% Inputs: 
% x      1D list of values you want on the x axis
% y      1D list of values you want on the y axis
% Z      2D grid of values you want to plot the feasible regions on
% c0     Inequality constraint Z > c0 or Z < c0
% c1     Displays a dotten contour next to the c0 contour - this signifies which side of the c0 contour you want to look. c1 should satisfy the Z constraint.
%
% No outputs
%
% V2.3.1, Eric Comstock, 1700 ET 02/04/2024 AD

% Inputs are the x, y (note these two should be 1D), and 2D Z values for the contour plot, the equality Z = c0, and a c1 on the right side of Z = c0 to show
% where the inequality should be. color is the color in a string
    cc0 = contourc(x, y, Z,[c0, c0]);                                                   % generates the points corresponding to Z = c0 - the region boundary for the inequality
    cc1 = contourc(x, y, Z,[c1, c1]);                                                   % generates the points corresponding to Z = c1 - showing what is inside the allowable region
    hold on                                                                             % allows multiple plots to be plotted at once
    plot(cc0(1,2:end), cc0(2,2:end), 'Color', color, 'LineStyle', '-', 'LineWidth', 3)  % plots the Z = c0 line with a thick solid line
    plot(cc1(1,2:end), cc1(2,2:end), 'Color', color, 'LineStyle', '-.', 'LineWidth', 3) % plots the Z = c1 line with a thick, dasked line for contrast
end

