function [r_mag,k,psi,Br_orient_strArr] = RotatingHB_magnetPose(a,b,n_mag,half_range)
%RotatingHB_magnetPose
%   Input:     a - parabolic curvature in x dimension
%                  b - parabolic curvature in y dimension  
%                  n_mag - number of magnets in one each dimension
%                  half_range - The last magnets location in half range
%
%  Output:   r_mag - 52 by 3 array, each row containst the position of one magnet
%                  k - single-rotation axis in cartesian coordinates
%                  psi - single-rotation angle [deg]
%                  Br_orient_strArr - Magnetization orientation in material frame [deg]


% Generate magnets position according to Halbach Array geometry
x = linspace(-half_range,half_range,n_mag);
y = linspace(-half_range,half_range,n_mag);
[X,Y] = meshgrid(x,y);
Z       = a*X.^2+b*Y.^2;

% Initialize the output array
r_mag     = zeros(n_mag*n_mag,3);
k_hat       = zeros(n_mag*n_mag,3);
Br_orient_strArr = zeros(n_mag*n_mag,3);
k = zeros(n_mag*n_mag,3);
psi = zeros(n_mag*n_mag,1);

for i = 1:(n_mag*n_mag)
    % Position of magnets
    r_mag(i,:)      = [X(i) Y(i) Z(i)];
    % Normal vector of the parabolic surface
    gradF = [-2*a*X(i),-2*b*Y(i),1];
    k_hat(i,:) = gradF/norm(gradF);
    % Find two sequential rotation angles required to align the z-axis of the magnet to the normal vector of the parabolic surface defining the parabolic shape of the halbach array 
    theta = atan2(-k_hat(i,2),k_hat(i,3));         % rotation angle around x-axis
    phi = atan2(k_hat(i,1),k_hat(i,3));              % rotation angle around y -axis

    % Combine two sequential rotations into one single rotation using Rodrigues' rotation formula
    Rx = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
    Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];
    R = Ry*Rx;
    
    % Find the single-rotation angle and single-rotation axis
    psi(i) = acos((trace(R)-1)/2);
    k(i,:) = [(R(3,2)-R(2,3))/2/sin(psi(i)),(R(1,3)-R(3,1))/2/sin(psi(i)),(R(2,1)-R(1,2))/2/sin(psi(i))];
    
    % Define the remanent flux density direction (Magnetization direction) in material frame
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

% Convert angle to degree
psi = rad2deg(psi);
end

