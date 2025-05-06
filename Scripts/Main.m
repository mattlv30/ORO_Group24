%% ORO Project Group 24

clear all; clc; close all; 

%% Assumptions 

%Earth radius and gravitational parameter
R_Earth = 6378.1370; %[km]
mu_Earth = 3.98600418e5; %[km^3/s^2]

% Altitude and radius vector of the spacecraft S/C
h1 = 600; %[km]
r1 = R_Earth + h1; %[km]

% Altitude and radius vector of the Debris 
h2 = 350; %[km]
r2 = R_Earth + h2; %[km]

%Phasing angle 
phi0_deg = 10; %[°]
phi0_rad = deg2rad(phi0_deg); %[rad]

% Orbital Velocity 
V1 = sqrt(mu_Earth/r1); %[km/s]
V2 = sqrt(mu_Earth/r2); %[km/s]

%% Hohmann Transfert - Analitical 

[DeltaV1,DeltaV2, TOF] = HohmannTransfer(r1,r2, mu_Earth);



