% Template for Assignment 3, part 2 in 30554
% (c) Daniel Olesen, DTU Space, 2023

% close all
clear
clc

addpath(genpath('GNSS-toolbox-master'))
addpath('navfunctions'); % Folder containing satpos + otherfunctions

% Input Rinex files. Remember that the observation file should be a Rinex 3.X format 
% and nav-file should be in version 2.XX  

% obs1 = 'static1.22O';
% nav1 = 'static1.22N';%
% obs1 = 'good\ALT30600.24O';
% nav1 = 'good\ALT30600.24N';
obs1 = 'kin\ALT30600.24O';
nav1 = 'kin\ALT30600.24N';

param = OBSRNX.getDefaults();
param.filtergnss = 'G';

% First time run (load from RINEX)
data1 = OBSRNX(obs1,param);

% Constants

f1=1575.42e6;
f2=1227.60e6;
c = 299792458; % m/s

% Import GPS ephemeris
brdc1 = loadRINEXNavigation('G','',nav1);
[eph,ion] = brdc2eph(brdc1,data1.t(500,8));

alfa = [ion.alpha0 ion.alpha1 ion.alpha2 ion.alpha3];  % From importGPSephemeris
beta = [ion.beta0 ion.beta1 ion.beta2 ion.beta3];      % From importGPSephemeris
%% Calculate GPS positions

% Processing options

elv_mask = 5; % Define which elevation mask to use (0 includes all SVs)       
useSaastamoinen = false; % [true | false];
useKlobuchar = false; % [true | false];
useIonofreeCombination = false; % [false | true]; Disable klobuchar if used

% Position calculations

c = 299792458; % m/s
f1=1575.42e6;
f2=1227.60e6;

x0=data1.recpos; %Use receiver position from Rinex header as 'truth'
lla = ecef2lla(x0);

lambda = lla(1);
phi = lla(2);
height = lla(3);

SV = data1.sat.G; 
pseudoranges_L1 = zeros(length(data1.t),length(SV));
pseudoranges_L2 = zeros(length(data1.t),length(SV));

for i=1:length(SV)
    pseudoranges_L1(:,i)=data1.obs.G{i}(:,1);
    pseudoranges_L2(:,i)=data1.obs.G{i}(:,5);
end

% Provide a function pr_ion_free which calculates the ionosphere free
% combination based on pseudorange measurements from L1 and L2 

% Select to use regular pseudoranges or ionospheric_free combination
if(useIonofreeCombination)
    pseudoranges_ionofree_L1 = pr_ion_free(pseudoranges_L1,pseudoranges_L2, f1,f2);
    pseudoranges = pseudoranges_ionofree_L1;
else
    pseudoranges = pseudoranges_L1;
end

prnList = SV;

% % exclude PRNs not in ephemeris
prnList([9])=[];
pseudoranges(:,[9])=[];

% define matrix for azimuth and elevatons
azimuth=zeros(size(pseudoranges));
elevation=zeros(size(pseudoranges));
pos = zeros(size(pseudoranges,2),3);
rec_d = zeros(size(pseudoranges,2),1);

% Calculate rotation-matrix for ENU conversion

R_L = R1(90-lambda)*R3(phi+90); % For conversion to ENU

lambda_L1 = 0.19029367;  % Wavelength for GPS L1 frequency in meters
lambda_L2 = 0.244210213; % Wavelength for GPS L2 frequency in meters

estimated_multipath = zeros(9, length(data1.t));

for i = 1:length(data1.t) % time-counter
    
    % find "active" observations  
    ix = (pseudoranges(i,:)>2e7 & pseudoranges(i,:)<3e7);
    pr = pseudoranges(i,ix);
    prn = prnList(ix);
    
    satPositions=zeros(3, length(pr));
    satClkCorr=zeros(1, length(pr));
    
    if length(pr) < 3
        continue
    end

    % Calculate satellite positions and account for earth rotation
        for j=1:length(satPositions)
    
                [satPositions(:,j), ~, satClkCorr(j)] = satposition(data1.t(i,8)-(pr(j)/c), prn(j), eph,useIonofreeCombination);
                % Account for earth-rotation
                satPositions(:,j) = e_r_corr(((pr(j)/c)+satClkCorr(j)), satPositions(:,j));    
        end
    
    % Find azimuth and elevaations for each satellite (Provide your own implementation here)   
    [azim, elev] = AzimuthElevation(R_L, satPositions, x0); 
    az(i,ix)=azim; el(i,ix)=elev;
    
    if(useSaastamoinen) % Executed if useSaastamoinen = true
        %Saastamoinen model
        for j=1:length(pr)
            dtrop = saastamoinen(elev(j), lambda, height/1000, 1013, 273+18, 0.5);
            pr(j) = pr(j) - dtrop;
        end
    end
    
    if(useKlobuchar) % Executed if useKlobuchar = true
       % Klobuchar model
        for j=1:length(pr)
            dIon = klobuchar(lambda,phi,elev(j),azimuth(j),data1.t(i),alfa,beta);
            pr(j) = pr(j) - dIon;
        end  
    end
    
    for j = 1:length(pr)
        p_L1 = data1.obs.G{j}(i, 1);
        phi_L1 = data1.obs.G{j}(i, 2) * lambda_L1
        phi_L2 = data1.obs.G{j}(i, 6) * lambda_L2
        if phi_L1 == 0 || phi_L2 == 0
            continue
        end
        MP_L1 = p_L1 - ((f1^2 + f2^2) / (f1^2 - f2^2)) * phi_L1 +  ((2 * f2^2) / (f1^2 - f2^2)) * phi_L2
        estimated_multipath(j, i) = MP_L1;
    end
    
    % Compensate for the unknown K1 by subtracting the mean of the estimated multipath
    mean_multipath = mean(estimated_multipath);
    corrected_multipath = estimated_multipath - mean_multipath;
    
    idx = elev > elv_mask;
    
    % [pos(i,:), rec_d(i)]= LeastSquaresGPS((pr(idx)+(c*satClkCorr(idx))), satPositions(:,idx));

end
%% PART 3.2
% Assuming `estimated_multipath` is already computed as described before

% Define colors for each satellite (optional)
colors = turbo(size(estimated_multipath, 1));

% Plot the estimated multipath effects for each satellite
figure;
hold on;
for j = 1:size(estimated_multipath, 1)
    plot(1:size(estimated_multipath, 2), estimated_multipath(j, :), 'Color', colors(j, :), 'LineWidth', 1.2);
end
hold off;

% Add labels and title
xlabel('Observation Index');
ylabel('Estimated Multipath (meters)');
title('Estimated Multipath Effects for Each Satellite');

% Add a legend (optional)
legendLabels = cell(length(prnList), 1);
for i = 1:length(prnList)
    legendLabels{i} = ['Satellite ', num2str(prnList(i))];
end
legend(legendLabels, 'Location', 'best');
grid on;

% You can customize this plot further by adding grid lines, labels, or any other decorations based on your preferences.
