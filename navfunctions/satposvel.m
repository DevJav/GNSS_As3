function [satPositions, satVelocity, satClkCorr] = satposvel(transmitTime, prnList, ...
                                             eph) 

%SATPOSVEL Computation of satellite coordinates X,Y,Z and Veocity at 
% TRANSMITTIME for ephemeris eph. Coordinates are computed for each satellite in the
%list PRNLIST.
%[satPositions, satClkCorr] = satpos(transmitTime, prnList, eph);
%
%   Inputs:
%       transmitTime  - transmission time
%       prnList       - list of PRN-s to be processed
%       eph           - ephemerides of satellites
%
%   Outputs:
%       satPositions  - position of satellites (in ECEF system [X; Y; Z;])
%       satClkCorr    - correction of satellite clocks
%
% Based on Kai Borre 04-09-96
% Copyright (c) by Kai Borre
% Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% Modified and extended to calculate satellite velocities by Daniel Olesen, DTU Space 2016
%
% Satellite velocity calculations is based on: http://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c
%
%% Initialize constants ===================================================
numOfSatellites = size(prnList, 2);

% GPS constatns

gpsPi          = 3.1415926535898;  % Pi used in the GPS coordinate 
                                   % system

%--- Constants for satellite position calculation -------------------------
Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM             = 3.986005e14;      % Universal gravitational constant times
                                   % the mass of the Earth, [m^3/s^2]
F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);
satVelocity = zeros(3, numOfSatellites);
sat_td = zeros(1, numOfSatellites);

%% Process each satellite =================================================

for satNr = 1 : numOfSatellites
    
    prn = prnList(satNr);
    prn = sprintf('PRN%02d', prn);
    
%% Find initial satellite clock correction --------------------------------

    %--- Find time difference ---------------------------------------------
    dt = check_t(transmitTime - eph.(prn).t_oc);
    
    
    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph.(prn).a_f2 * dt + eph.(prn).a_f1) * dt + ...
                         eph.(prn).a_f0 - ...
                         eph.(prn).T_GD;

    time = transmitTime - satClkCorr(satNr);

%% Find satellite's position ----------------------------------------------

    %Restore semi-major axis
    a   = eph.(prn).sqrtA * eph.(prn).sqrtA;

    % 1) Time correction
    tk  = check_t(time - eph.(prn).t_oe);

    %Initial mean motion
    n0  = sqrt(GM / a^3);
    %Mean motion
    n   = n0 + eph.(prn).deltan;

    % 2) Mean anomaly
    M   = eph.(prn).M_0 + n * tk;
    
    % Derivative of Mean anomaly
    
    MDot = n;
    
    %Reduce mean anomaly to between 0 and 360 deg
    
    M   = rem(M + 2*gpsPi, 2*gpsPi);

    %Initial guess of eccentric anomaly
    E   = M;

    % 3) Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph.(prn).e * sin(E);
        
       % EDot = MDot / (1.0 - eph.(prn).e*cos(E));
        
        dE      = rem(E - E_old, 2*gpsPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end

    EDot = MDot / (1.0 - eph.(prn).e*cos(E));
       
    
    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*gpsPi, 2*gpsPi);

    %Compute relativistic correction term
    dtr = F * eph.(prn).e * eph.(prn).sqrtA * sin(E);

    % 4) Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph.(prn).e^2) * sin(E), cos(E)-eph.(prn).e);

    %Calculate the derivative of the true anomaly
    nuDot = sin(E)*EDot*(1.0 + eph.(prn).e*cos(nu))/( sin(nu)*(1-eph.(prn).e*cos(E)));
    
    %5) Compute angle phi (latitude)
    phi = nu + eph.(prn).omega;
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*gpsPi);
    
    corr_u = eph.(prn).C_us*sin(2.0*phi) + eph.(prn).C_uc*cos(2.0*phi);
    corr_r = eph.(prn).C_rs*sin(2.0*phi) + eph.(prn).C_rc*cos(2.0*phi);
    corr_i = eph.(prn).C_is*sin(2.0*phi) + eph.(prn).C_ic*cos(2.0*phi);
    
    u = phi + corr_u;
    r = a * (1 - eph.(prn).e*cos(E)) + corr_r;
    i = eph.(prn).i_0 + eph.(prn).iDot * tk + corr_i;
    

%     %Correct argument of latitude
%     u = phi + ...
%         eph.(prn).C_uc * cos(2*phi) + ...
%         eph.(prn).C_us * sin(2*phi);
%     %Correct radius
%     r = a * (1 - eph.(prn).e*cos(E)) + ...
%         eph.(prn).C_rc * cos(2*phi) + ...
%         eph.(prn).C_rs * sin(2*phi);
%     %Correct inclination
%     i = eph.(prn).i_0 + eph.(prn).iDot * tk + ...
%         eph.(prn).C_ic * cos(2*phi) + ...
%         eph.(prn).C_is * sin(2*phi);
    
    % Calculate derivative of latitude, radius and inclination
    
    uDot = nuDot + 2*(eph.(prn).C_us*cos(2*u)-eph.(prn).C_uc*sin(2*u))*nuDot;
    rDot = (a*eph.(prn).e*sin(E)*n)/(1-eph.(prn).e*cos(E)) + ...
        2*(eph.(prn).C_rs*cos(2*u)-eph.(prn).C_rc*sin(2*u))*nuDot;
    iDot = eph.(prn).iDot + (eph.(prn).C_is*cos(2*u)-eph.(prn).C_ic*sin(2*u))*2*nuDot;
    
    

    %Compute the angle between the ascending node and the Greenwich meridian
  Omega = eph.(prn).omega_0 + (eph.(prn).omegaDot - Omegae_dot)*tk - ...
            Omegae_dot * eph.(prn).t_oe;

    %Reduce to between 0 and 360 deg
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);
    
    % Compute the derivate of The Angle
    
    OmegaDot = (eph.(prn).omegaDot - Omegae_dot);
    
    % Define intermediate calculations for simplicity
    
    xpk = r*cos(u); 
    ypk = r*sin(u); 
    xpkDot = rDot*cos(u) - ypk*uDot;
    ypkDot = rDot*sin(u) + xpk*uDot;

    %--- Compute satellite coordinates ------------------------------------
    satPositions(1, satNr) = cos(u)*r * cos(Omega) - sin(u)*r * cos(i)*sin(Omega);
    satPositions(2, satNr) = cos(u)*r * sin(Omega) + sin(u)*r * cos(i)*cos(Omega);
    satPositions(3, satNr) = sin(u)*r * sin(i);
    
    %--- Compute satellite velocities ------------------------------------
    satVelocity(1, satNr) = ( xpkDot-ypk*cos(i)*OmegaDot )*cos(Omega)...
        - ( xpk*OmegaDot+ypkDot*cos(i)-ypk*sin(i)*iDot )*sin(Omega);
    
    satVelocity(2, satNr) = ( xpkDot-ypk*cos(i)*OmegaDot )*sin(Omega)...
        + ( xpk*OmegaDot+ypkDot*cos(i)-ypk*sin(i)*iDot )*cos(Omega);
    
    satVelocity(3, satNr) = ypkDot*sin(i) + ypk*cos(i)*iDot;


%% Include relativistic correction in clock correction --------------------
    satClkCorr(satNr) = (eph.(prn).a_f2 * dt + eph.(prn).a_f1) * dt + ...
                         eph.(prn).a_f0 - ...
                         eph.(prn).T_GD + dtr;
                                          
end % for satNr = 1 : numOfSatellites
