function [X,Y,F] = spacingconv(X,Y,x,y,f)
% Returns the input data into appropriate format for gridinterp_Pi_m
% [X,Y,F] = spacingconv(X,Y,x,y,f)
%
% Inputs:
%     1st Input  : X, ferrofluid x coordinate array
%     2nd Input  : Y, ferrofluid y coordinate array
%     3rd Input  : x, expanded X (column vector)
%     4th Input  : y, expanded Y (column vector)
%     5th Input  : f, interpolated data from COMSOL
%
% Outputs:
%     X: 1st Input
%     Y: 2nd Input
%     F: Rearranged f 
%
% V1.1 - Álvaro Romero-Calvo (04/01/2024)

F     = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        F(i,j) = f(x == X(i,j) & y == Y(i,j));
    end
end
end