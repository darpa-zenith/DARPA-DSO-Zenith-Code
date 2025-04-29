function [Fk_x,Fk_y,Fk] = Kelvin_H(M_vec, Hx,Hy,params)
    %% Kelvin_F calculates the Kelvin Body force from magnetic field H
    %
    % Inputs: 
    %       M_vec          
    %       Hx            
    %       Hy
    %       params
    %
    % Outputs: 
    %       Fk_x  Magnetic mass-specific Kelvin body force, x comp [m/s^2]
    %       Fk_y  Magnetic mass-specific Kelvin body force, y comp [m/s^2]
    %       Fk    Magnetic mass-specific Kelvin body force [m/s^2]
    %
    % V1.0.1, Hugh Chen, 05/04/2024 AD
    
    %% Parameters
    hx_eric = params.hx_eric;
    hy_eric = params.hy_eric;
    mu0     = params.mu0;
    rho     = params.rho;

    %% Calculate Kelvin Body Force
    [Hx_x,Hx_y]    = gradient(Hx,hx_eric,hy_eric);   %intermediate step for calculation of gradient
    [Hy_x,Hy_y]    = gradient(Hy,hx_eric,hy_eric);   %intermediate step for calculation of gradient
    Fk_x           = mu0 * (dot(M_vec,cat(3,Hx_x,Hx_y),3)); %Kelvin body force in x dir
    Fk_y           = mu0 * (dot(M_vec,cat(3,Hy_x,Hy_y),3)); %Kelvin body force in y dir  
    Fk             = sqrt(Fk_x.^2+Fk_y.^2)/rho/9.81;        %[m/s^2] Kelvin body force magnitude, in g
end