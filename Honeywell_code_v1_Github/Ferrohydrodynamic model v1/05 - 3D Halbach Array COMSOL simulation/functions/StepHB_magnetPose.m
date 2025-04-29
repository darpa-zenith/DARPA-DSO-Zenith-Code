function [r_mag,Br_orient_strArr] = StepHB_magnetPose(a,b,n_mag,half_range)
%StepHB_magnetPose
%   Input:     a - parabolic curvature in x dimension
%                  b - parabolic curvature in y dimension  
%                  n_mag - number of magnets in one each dimension
%                  half_range - The last magnets location in half range
%
%  Output:   r_mag - 52 by 3 array, each row containst the position of one magnet
%                  Br_orient_strArr - Magnetization orientation in material frame [deg]

% Generate magnets position according to Halbach Array geometry
x = linspace(-half_range,half_range,n_mag);
y = linspace(-half_range,half_range,n_mag);
[X,Y] = meshgrid(x,y);
Z       = a*X.^2+b*Y.^2;

% Initialize the output array
r_mag     = zeros(n_mag*n_mag,3);
Br_orient_strArr = zeros(n_mag*n_mag,3);

for i = 1:(n_mag*n_mag)
    % Position of magnets
    r_mag(i,:)      = [X(i) Y(i) Z(i)];

    % Define the remanent flux density direction (Magentization direction) in material frame
    if mod(i,4)==1
        Br_orient_strArr(i,:) = string([0 0 1]);
    elseif mod(i,4)==2
        Br_orient_strArr(i,:) = string([0 -1 0]);
    elseif mod(i,4)==3
        Br_orient_strArr(i,:) = string([0 0 -1]);
    elseif mod(i,4)==0
        Br_orient_strArr(i,:) = string([0 1 0]);
    end
end
end

