% Template for Assignment 3, part 1 in 30554
% (c) Daniel Olesen, DTU Space, 2023

close all
clear
clc

addpath(genpath('GNSS-toolbox-master'))
addpath('navfunctions'); % Folder containing satpos + otherfunctions

% Input Rinex files. Remember that the observation file should be a Rinex 3.X format 

obs1 = 'kin\ALT30600.24O';
nav1 = 'kin\ALT30600.24N';

param = OBSRNX.getDefaults();
param.filtergnss = 'G';

% load data from RINEX
data1 = OBSRNX(obs1,param);

% Plot measurements 

%% ----------  Fill in your own code for analyzing data -------------- %

for i=1:11
    obs = data1.obs.G(i);
    obs_data = obs{1};
    C1C = obs_data(:,1);
    L1C = obs_data(:,2);
    D1C = obs_data(:,3);
    S1C = obs_data(:,4);

    figure(1);
    hold on;
    plot(C1C, 'DisplayName', ['Satellite ', num2str(data1.sat.G(i))]);
    title("C1C");
    xlabel("Time");
    ylabel("Pseudoranges (m)");
    grid();

    figure(2);
    hold on;
    plot(L1C, 'DisplayName', ['Satellite ', num2str(data1.sat.G(i))]);
    title("L1C");
    xlabel("Time");
    ylabel("Carrier phase");
    grid();

    figure(3);
    hold on;
    plot(D1C, 'DisplayName', ['Satellite ', num2str(data1.sat.G(i))]);
    title("D1C");
    xlabel("Time");
    ylabel("Doppler (Hertz)");
    grid();

    figure(4);
    hold on;
    plot(S1C, '.', 'DisplayName', ['Satellite ', num2str(data1.sat.G(i))])
    title("S1C");
    xlabel("Time");
    ylabel("Carrier-to-noise (C/N0)");
    grid();
end

% Adding legends
for j = 1:4
    figure(j);
    legend('show');
end
