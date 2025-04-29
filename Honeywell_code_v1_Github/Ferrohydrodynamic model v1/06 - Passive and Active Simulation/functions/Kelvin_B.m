function [Fk] = Kelvin_B(x,z,params)
    %% Kelvin_F calculates the Kelvin Body force from magnetic flux density B 
    %
    % Inputs: 
    % x             simulation x coord                  
    % z:            simulation z/y coord                
    % params        Parameter variable structure
    %
    % Outputs: 
    % Fk            Magnetic mass-specific Kelvin body force             [m/s^2]
    %
    % V2.0.1, Hugh Chen, 21/04/2024 
    % V1.0.1, Eric Comstock, 01/04/2024 AD
    

    % Parameters
    eps_     = 1e-10;
    Pi_m_gi  = params.Pi_m;
    %% Calculating Mass Force Potential FDM
    x        = x';
    z        = z';
    dPi_m_dx = (Pi_m_gi(x+eps_,z)-Pi_m_gi(x-eps_,z))/2/eps_;    % x derivative of magnetic potential
    dPi_m_dz = (Pi_m_gi(x,z+eps_)-Pi_m_gi(x,z-eps_))/2/eps_;    % z derivative of magnetic potential

    %Converting back to original position
    dPi_m_dx = dPi_m_dx';       
    dPi_m_dz = dPi_m_dz';

    Fk_x     = -dPi_m_dx;
    Fk_z     = -dPi_m_dz;                   
    Fk       = sqrt(Fk_x.^2 + Fk_z.^2);  % KBF magnitude
end