function p_new = extractpoints(p,sphere_center,distance2SphereCenter)
    % Use an invisible scan sphere located at the focal point to extract the top parabolic curve
    % of the equipotential lines
    %
    % Input:    p - old data points of the contour line
    %               sphere-center - Scan Sphere center (focal point)
    %               distance2SphereCenter - invisible Scan Sphere radius [m]
    %
    %Ouput:    p_new - new data points of the contour line

    r = vecnorm(sphere_center-p);                                                   % Find distance from all data points to the focal point 
    p_new = p(:,r<=distance2SphereCenter & abs(p(1,:))<=0.25);   % Find the datapoints whose distance to the focal point is smaller that the input distance2SphereCenter
    [p_new(1,:),Ind] = sort(p_new(1,:),2);                                              % Sort the datapoints in ascending x order
    p_new(2,:) = p_new(2,Ind);
end
