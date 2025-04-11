function val = interpXYZ(Xin, Yin, Zin, Hin, Xout, Yout, Zout)
%INTERPXYZ
%   Xin, Yin, Hin,  input Hin(Xin, Yin, Zin)
%   return pts = Hin(Xout, Yout, Zout)

% This assumes Xin and Yin values are sorted in ascending order and that 
% Xin(j,:) contains a single value and Yin(:,k) all contains the same value

numx = size(Hin, 1);
numy = size(Hin, 2);
numz = size(Hin, 3);

if numx ~= size(Xin, 1)
    disp("X sizes vary");
end

if numy ~= size(Yin, 2)
    disp("Y sizes vary");
end

if numz ~= size(Zin, 3)
    disp("Z sizes vary");
end

% assume angles are same for all x and z
AngOut = atan2(Yout(1,:), Xout(1,:));
if AngOut < 0
    AngOut = AngOut + 2*pi;
end

% First interpolate in r and theta directions simultaneously, in two different
% planes.  Then interpolate between the points in those two z-planes with
% constant (r, theta) coordinates.

% Assume z-values are constant across x and y
% Same for x and y

xxs = zeros(numx, 1);
yys = zeros(numy, 1);
zzs = zeros(numz, 1);

% radial dimension
for i = 1:size(Xin,1)
    xxs(i) = norm([Xin(i,1,1), Yin(i,1,1)]);
end

% angle (radians)
for j = 1:size(Yin,2)
    yys(j) = atan2(Yin(1,j,1), Xin(1,j,1));
    if yys(j) < 0
        yys(j) = yys(j) + 2*pi;
    end
end

% height
for i = 1:size(Zin,3)
    zzs(i) = Zin(1,1,i);
end

val = zeros(length(Xout), 1);
for ii=1:length(Xout)

% ----------------------------------------    

    indz1 = find(zzs > Zout(ii), 1);
    indz0 = indz1 - 1;

    if zzs(indz1) - zzs(indz0) < 0.0000001
        omega = 1;
    else
        omega = (Zout(ii) - zzs(indz0)) / (zzs(indz1) - zzs(indz0));
    end

% ----------------------------------------   

    r = norm([Xout(ii), Yout(ii)]);
    indr = find(xxs > r, 1);

    if numy == 1

        omega1 = (r - xxs(indr-1)) / (xxs(indr) - xxs(indr-1));
        val0 = (1-omega1)*Hin(indr-1,1,indz0) + omega1*Hin(indr,1,indz0);
        val1 = (1-omega1)*Hin(indr-1,1,indz1) + omega1*Hin(indr,1,indz1);

    else

        angle = atan2(Yout(ii), Xout(ii));
        if angle < 0
            angle = angle + 2*pi;
        end
        indtheta = find(yys > angle, 1);

        denom = 1.0/((xxs(indr)-xxs(indr-1)*(yys(indtheta)-yys(indtheta-1))));
        xvec = [xxs(indr,1)-r  r-xxs(indr-1,1)];
        yvec = [yys(indtheta)-angle; angle-yys(indtheta-1)];
    
        matxy0 = [ Hin(indr-1, indtheta-1, indz0)   Hin(indr-1, indtheta, indz0);
                   Hin(indr, indtheta-1, indz0)     Hin(indr, indtheta, indz0)] ;
    
        val0 = denom * xvec * matxy0 * yvec;

        matxy1 = [ Hin(indr-1, indtheta-1, indz1)   Hin(indr-1, indtheta, indz1);
                   Hin(indr, indtheta-1, indz1)     Hin(indr, indtheta, indz1)] ;
        val1 = denom * xvec * matxy1 * yvec;
    end

% ----------------------------------------  

    val(ii,1) = (1-omega)*val0 + omega*val1;

% ----------------------------------------    

end

end