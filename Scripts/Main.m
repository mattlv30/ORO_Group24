%% ORO Project Group 24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; 

%% Initial data 

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
n1=sqrt(mu_Earth/r1^3); 
n2=sqrt(mu_Earth/r2^3); 

%Orbital Element
e = 0;
i = 0;
w = 0;
om =0;

%% Part 1 - Hohmann transfer

[DeltaV1,DeltaV2, ToF, phi, t_phasing] = part1(r1,r2, mu_Earth, phi0_rad, V1, e , i, w, om);

%% part 2
[rH_err,rH_dot_err,Delta_T,deltaV1,deltaV2] = part2(r2,n2);

aaa