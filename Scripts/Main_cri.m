%% ORO Project Group 24

clear; clc; close all; 

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

% Orbital Velocity 
V1 = sqrt(mu_Earth/r1); %[km/s]
V2 = sqrt(mu_Earth/r2); %[km/s]
n2=sqrt(mu_Earth/r2^3);

%% part 2
[rH_err,rH_dot_err,Delta_T,deltaV1,deltaV2] = part2(r2,n2);

%% part 3

rho=7.014*10^(-12)*10^9; % atmospheric density [kg/m^3] ---> [kg/km^3]
CD=2.2; % drag coefficient
A=2/10^6; % cross-sectional area sat [m^2] ---> [km^2]
m=325; % sat mass [kg]

part3(rho,CD,A,m,Delta_T,r2,mu_Earth,rH_err,rH_dot_err,deltaV1,deltaV2)

%% part 5

pre_part5(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)

part5_vers1(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)

part5_vers2(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)