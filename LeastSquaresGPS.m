function [pos,clock_err] = LeastSquaresGPS(R, satellites)
    
    % initiate preliminary position
    x_0 = [0,0,0,0];
    c = 299792458; % m/s
    num_satellites = length(satellites);

    distance_to_satellites = zeros(num_satellites, 1);
    A = zeros(num_satellites, 4);
    z = zeros(num_satellites, 1);

    delta_x = [1000,1000,1000,1000];
    iter = 0;

    while abs(delta_x(1)) > 1e-10 || iter > 10
        for i=1:num_satellites
            sat = satellites(:,i);
            distance_to_satellites(i) = sqrt((sat(1) - x_0(1))^2 + (sat(2) - x_0(2))^2 + (sat(3) - x_0(3))^2);

            A(i, 1) = - ((sat(1) - x_0(1)) / distance_to_satellites(i));
            A(i, 2) = - ((sat(2) - x_0(2)) / distance_to_satellites(i));
            A(i, 3) = - ((sat(3) - x_0(3)) / distance_to_satellites(i));
            A(i, 4) = 1;

            z(i) = R(i) - distance_to_satellites(i) - c * x_0(4);
        end
        delta_x = A\z;
        
        x_0 = x_0 + delta_x;
        iter = iter + 1;
    end

    pos = x_0(1:3);
    clock_err = x_0(4);

end