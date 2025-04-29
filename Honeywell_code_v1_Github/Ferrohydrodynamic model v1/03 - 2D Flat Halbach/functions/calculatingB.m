function [Bx,By] = calculatingB(params,h,w,X,Y,theta,r0)
% This function calculates the two magnetic flux density vector components
% [Bx, By] in 2D. 
%
%  Input  : params
%           h     (Heights of magnets)                  [m]
%           w     (Widths of magnets)                   [m]
%           X     (x coordinates matrix of Mesh points) [m]
%           Y     (y coordinates matrix of Mesh points) [m]
%           theta (Magnets orientations with respect to
%                  vertical axis)                       [rad] 
%           r0    (Position vectors of each magnets)    [m]
%
%  Output : Bx   (x components of magnetic flux density vectors) [T]
%           By   (y components of magnetic flux density vectors) [T]



    %% Parameters
    mu0 = params.mu0;
    Ke = params.Ke;

    %% Calculation
    Bx        = zeros(size(X));
    By        = zeros(size(X));

    for i = 1:size(r0,1)
        xi    = cos(theta(i)) * (X-r0(i,1)) + sin(theta(i)) * (Y - r0(i,2));
        psi   = -sin(theta(i)) * (X-r0(i,1)) + cos(theta(i)) * (Y - r0(i,2));

        B_xi  = mu0*Ke(i)/2/pi/2*(-log(((xi-w(i)/2).^2+(psi-h(i)/2).^2)./((xi-w(i)/2).^2+(psi+h(i)/2).^2))+log(((xi+w(i)/2).^2+(psi-h(i)/2).^2)./((xi+w(i)/2).^2+(psi+h(i)/2).^2)));
        B_psi = mu0*Ke(i)/2/pi*(-(atan((psi+h(i)/2)./(xi-w(i)/2)))+atan((psi-h(i)/2)./(xi-w(i)/2))+(atan((psi+h(i)/2)./(xi+w(i)/2)))-atan((psi-h(i)/2)./(xi+w(i)/2)));

        Bx    = Bx+cos(theta(i))*B_xi-sin(theta(i))*B_psi;
        By    = By+sin(theta(i))*B_xi+cos(theta(i))*B_psi;
    end
end











