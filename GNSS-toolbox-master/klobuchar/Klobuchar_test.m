clear all
clc
%change the path
addpath(genpath('C:\Users\saschu\Danmarks Tekniske Universitet\Geo_GNSS - Dokumenter\Teaching\30554 GNSS (new course)\GNSS-toolbox-master\src'));

param = OBSRNX.getDefaults();
param.filtergnss = 'G';

%change the file names if needed
obsfileGRL = 'qaq2307j.21O';
navfileGRL = 'qaq2307j.21N';
obsfileDK = '20220214_b328_roof.22O';
navfileDK = '20220214_b328_roof.22N';
obs = OBSRNX(obsfileGRL,param);

%% Load nav file & Klobuchar input
nav = loadRINEXNavigation('G','',navfileGRL);
% Input data needed for Klobuchar model
fi = -46.0479;             % Receiver's latitude [deg] the current latitude is for the station qaq2
lambda = 60.7153;         % Receiver's longitude [deg] the current longitude is for the station qaq2
elev = ;                     % Elevation angle of satellite [deg] ... as calculated in Assignment 1
azimuth = ;                   % Geodetic azimuth of satellite [deg]... as calculated in Assignment 1              
tow = obs.t(1,8);                     % Extract initial Time of Week 
alpha = nav.hdr.ionoAlpha;
beta = nav.hdr.ionoBeta;


%% Implement Klobuchar model --- Assignment 3 Part 2 
dIon = klobuchar(fi,lambda,elev,azimuth,tow,alpha,beta)% Ionospheric range correction for the GPS L1 frequency

