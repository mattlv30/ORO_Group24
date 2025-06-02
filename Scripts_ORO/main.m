%% ORO Project Group 24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all; 

folder = fullfile(pwd);
addpath(genpath(folder))

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

%% Part 2 - Relative dynamics final approach

[rH_err,rH_dot_err,Delta_T,deltaV1,deltaV2] = part2(r2,n2);

%% Part 3 - Perturbations 

rho=7.014*10^(-12)*10^9; % atmospheric density [kg/m^3] ---> [kg/km^3]
CD=2.2; % drag coefficient
A=2/10^6; % cross-sectional area sat [m^2] ---> [km^2]
m=325; % sat mass [kg]

 
part3(rho,CD,A,m,Delta_T,r2,mu_Earth,rH_err,rH_dot_err,deltaV1,deltaV2)

%% Part 4 - Robotic manipulator

% Contact force on the end effector 
F = 20; %[N]

% Wrench on the fifth link
w_5 = [0 0 0 0 -F 0]';


% Link geometry
L1 = 0.05; %[m]
L2 = 0.30; %[m]
L3 = 0.05; %[m]
L4 = 0.30; %[m]
L5 = 0.30; %[m]

L = [L1 L2 L3 L4 L5];

% Joint rotation

% 1, 3 and 5 are fixed in order to ensure the 2D motion of the manipulator 
theta1_deg = 0; %[°] 
theta1 = deg2rad(theta1_deg); %[rad]
theta3_deg = 0; %[°] 
theta3 = deg2rad(theta3_deg); %[rad]
theta5_deg = 0; %[°] 
theta5 = deg2rad(theta5_deg); %[rad]

% 2 and 4 can vary to consider different configurations
theta2_deg = 60; %[°] 
theta2 = deg2rad(theta2_deg); %[rad]
theta4_deg = 60; %[°] 
theta4 = deg2rad(theta4_deg); %[rad]

theta = [theta1 theta2 theta3 theta4 theta5];

% Solution of part 4

[tau] = part4(L,theta, w_5);




%% Part 5 - Optimisation of Delta V and ToF

pre_part5(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)

part5_vers1(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)

deltaVtot_vers2=part5_vers2(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1);

part5_vers3(rho,CD,A,m,r2,n2,rH_err,rH_dot_err,deltaV1,deltaVtot_vers2)