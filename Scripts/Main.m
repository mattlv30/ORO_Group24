clear all; clc; close all; 

%% Assumptions 

%Earth radius and gravitational parameter
R_Earth = 6378.1370; %[km]
mu_Earth = 3.98600418e5; %[km^3/s^2]

% Altitude and radius vector of the spacecraft S/C
h1 = 600; %[km]
R1 = R_Earth + h1; %[km]

% Altitude and radius vector of the Debris 
h2 = 350; %[km]
R2 = R_Earth + h2; %[km]

%Phasing angle 
phi0_deg = 10; %[°]
phi0_rad = deg2rad(phi0_deg); %[rad]

% Orbital Velocity 

V1 = sqrt(mu_Earth/R1); %[km/s]
V2 = sqrt(mu_Earth/R2); %[km/s]

%% Hohmann Transfert 

DeltaV1 = sqrt(mu_Earth/R1)*(sqrt(2*R2/(R1+R2))-1);
DeltaV2 = sqrt(mu_Earth/R2)*(1-sqrt(2*R2/(R1+R2)));



