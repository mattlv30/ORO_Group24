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

[DeltaV1,DeltaV2, DeltaVtot, a_h, e_h, ToF] = HohmannTransfer(r1,r2, mu_Earth);

fprintf('DeltaV1: %.4f km/s\n', DeltaV1);
fprintf('DeltaV2: %.4f km/s\n', DeltaV2);
fprintf('Total DeltaV: %.4f km/s\n', DeltaVtot);
fprintf('Time of Flight: %.1f s (%.2f minutes)\n', ToF, ToF/60);

%% Plotting the Orbits

span = 500;

% Circular orbits 
f = linspace(0, 2*pi, span); 
orbit1 = r1 * [cos(f); sin(f)];
orbit2 = r2 * [cos(f); sin(f)];

% Transfer ellipse
f_h = linspace(pi, 2*pi, span/2);
r_trans = a_h*(1-e_h^2) ./ (1+e_h*cos(f_h));
Htransfert = r_trans .* [cos(f_h) ; sin(f_h)];

figure;
hold on; grid on; axis equal;
plot(orbit1(1,:), orbit1(2,:), 'b--', 'LineWidth', 1.5);
plot(orbit2(1,:), orbit2(2,:), 'r--', 'LineWidth', 1.5);
plot(Htransfert(1,:), Htransfert(2,:), 'k-', 'LineWidth', 2);
plot(-r1, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b');
plot(r2*cos(2*pi), r2*sin(2*pi), 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
legend('Initial Orbit', 'Target Orbit', 'Hohmann Transfer', 'Start Point', 'End Point');
title('Hohmann Transfer Maneuver');
xlabel('X [km]'); ylabel('Y [km]');


%% fasatura
n2=sqrt(mu_Earth/r2^3);
n1=sqrt(mu_Earth/r1^3);
t_phasing=(n2*ToF-pi-phi0_rad)/(n1-n2)
phi=phi0_rad-(n2-n1)*t_phasing;
ToF*n2-phi