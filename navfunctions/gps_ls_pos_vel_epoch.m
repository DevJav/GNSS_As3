function [ pos, ctr, vel, dop ] = gps_ls_pos_vel_epoch( pos_init, Sat_Pos, Sat_Vel, Pseudo_range, Pseudo_range_rate, min_elevation, tow, ion, satClkCorr )
%GPS_LS calculates the gps positions from pseudorange and pseudo-range rate
% measurements using a LS approach
% 
% Inputs:
%   pos_init [3x1]: is the initial position vector(start-guess) 
%       given in curvilinear coordinates (lat, long, height)
%
%   sat_pos :                   Satellite positions (one row for each
%                               satellite)
%   sat_vel:                    Satellite velocity
%   pseudo_range:               Pseudoranges
%   pseudo_range_rate           Pseudorange-rates
%   
%   min_elevation  :            minimum elevation to include in computation
%
% Outputs:
%   pos                         ECEF user position
%   vel                         ECEF user velocity
%   DOP                         Dillution-of-Precision

% Author:
%    Daniel Olesen, DTU Space 2016
% Rev date: 2016-04-08

% Convert (initial) curvilinear coordinates to ECEF 

c = 299792458; % m/s

[x_i, y_i, z_i] = geodetic2ecef(pos_init(1)*(pi/180),pos_init(2)*(pi/180)...
    ,pos_init(3),wgs84Ellipsoid);

x_nom = [x_i y_i z_i 0 0 0 0 0]'; % Initial position, clock-bias, velocity and clock-drift
%x_nom = [x_i y_i z_i 0]';

% find all pseudo range and range rate measurements for first epoch

x = x_nom;

iter = 0; itermax = 20;
tol = 1e-12;
    
h = realmax*ones(size(x));
    
while (max(abs(h(1:4))) > tol) && (iter < itermax)

% Calculate satellite elevation angles
    
    noSatellites = length(Pseudo_range);
    
    w = zeros(1,noSatellites); % weighting vector for observations
    
    idx2 = []; % index for satellites below elevation mask

    [phi, lambda, h] = ecef2geodetic(x(1), x(2), x(3), wgs84Ellipsoid,'radians');    
    
    pr = Pseudo_range;
    prr = Pseudo_range_rate;
    X = Sat_Pos;
    VX = Sat_Vel;

    for i = 1:noSatellites
  
    % calculate approx travel-time 

    rho = sqrt((Sat_Pos(i,1)-x(1))^2 + (Sat_Pos(i,2)-x(2))^2 + (Sat_Pos(i,3)-x(3))^2); % geometric model
    traveltime = rho / c;

    %traveltime = pr(i) / c;

    % Calculate Updated satellite positions (at emission time)

    X(i,:) = e_r_corr(traveltime+satClkCorr(i), Sat_Pos(i,:)')';
    VX(i,:) = e_r_corr(traveltime+satClkCorr(i), Sat_Vel(i,:)')';

    
%     [xEast, yNorth, zUp] = ecef2enu(X(i,1), X(i,2), ...
%         X(i,3), phi, lambda, h, wgs84Ellipsoid,'radians'); 
%     
%     zenith = (360/(2*pi))* acos(zUp /...
%         ( sqrt( xEast^2 + yNorth^2 +  zUp^2 ) ) );
%        
%     %zUp / ( sqrt( xEast^2 + yNorth^2 +  zUp^2 )) 
%     
%     elevation_angle = 90-zenith
%     
%     azimuth = 180 - atan(xEast/yNorth)*(180/pi)
%   
    [az, el, dist] = topocent(x(1:3), X(i,:)' - x(1:3));


    
    if (~isempty(ion)) % if ionospheric parameters are given calculate correction according to the klobuchar model.
        [ ion_delay ] = iono_corr( el*(pi/180), az*(pi/180), phi, lambda, h, ion, tow );
        %ion_delay  = 0;   
    else
       ion_delay = 0; 
    end
    
    [dtrop] = saastamoinen(el, phi*(180/pi), h, 1013, 293, 0.5);
    %dtrop = 0;
    pr(i) = pr(i) - (ion_delay*c) - dtrop + satClkCorr(i)*c;

%     w(i) = 1;
     w(i) = (1 /(sind(el)^2));

       if(el < min_elevation)
           idx2 = [idx2; i]; % delete index if elevation is too low
       end
     end
    
    pr(idx2) = [];             % remove satellites with low elevations. 
    prr(idx2) = [];
    X(idx2,:) = []; VX(idx2,:) = [];
    w(idx2) = [];

   
    % Least-Squares (Gauss Newton) See Assignment F (Satellite-based
    % positioning course)
    
                
            Sat_distance = sqrt(sum((( X - repmat(x(1:3)',length(pr),1) ).^2)')');
            a_vec = (X - repmat(x(1:3)', length(pr),1))./ repmat(Sat_distance,1,3);
            
            z_k = [pr; prr];
            %z_k = [pr];

            z_k(1:length(pr)) = z_k(1:length(pr)) -  Sat_distance - (repmat(x(4)', length(pr),1));
            z_k(length(pr)+1:end) = z_k(length(pr)+1:end) - sum((VX.*a_vec)')';
    
%           z_k = z_k -  Sat_distance - (repmat(x(4)', length(pr),1));

            
            H_k = [-a_vec repmat(1,length(pr),1)];
             
            H_k = [H_k zeros(size(H_k));
                   zeros(size(H_k)) H_k];

            J = H_k;
         
            h= (sqrt(diag(1./[w w]))*J)\(sqrt(diag(1./[w w]))*z_k);
          
            %h= ((diag(1./w))*J)\((diag(1./w))*z_k);

            %h= (sqrt(diag(1./w))*J)\(sqrt(diag(1./w))*z_k);
                
            %h= J\z_k;   

            x = x + [h(1:4); zeros(4,1)];

            iter = iter + 1;
        end

        pos = x(1:3);
        ctr = x(4);
        vel = h(5:7);
        dop = 0;


    end

%end

