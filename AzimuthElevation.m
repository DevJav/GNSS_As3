function [az, el] = AzimuthElevation(R_L, satPos, x0)                                     
                                         
    for j=1:size(satPos,2)

    % Calculate ENU coordinates of SV (to calculate azimuth and Zenith)
    x_ENU = R_L*[satPos(1,j) - x0(1); satPos(2,j) - x0(2);satPos(3,j) - x0(3)];

    % Calculate Azimuth and Zenith
    azimuth = (360/(2*pi))*atan2(x_ENU(1),x_ENU(2));
    zenith = (360/(2*pi))*acos(x_ENU(3)/(sqrt(x_ENU(1)^2 + x_ENU(2)^2 + x_ENU(3)^2)));    

            if azimuth < 0
                azimuth = 360 + azimuth;
            end

            if ((90 - zenith) > 0)
                az(j)=azimuth;
                el(j)=90-zenith; % 
            else
                 az(j)=NaN;
                 el(j)=NaN; 
            end
    end

end