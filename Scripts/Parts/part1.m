%% Part 1
% This function solves all the request of the first part of the project

function [DeltaV1,DeltaV2, ToF, phi, t_phasing] = part1(r1,r2, mu_Earth, phi0_rad, V1, e , i, w, om)

%% Inputs
% r1 - raadius of the chaser orbit
% r2 - radius of the target orbit
% mu_Earth - gravitational pararmeter of the primary (Earth)
% phi_zero - initial phasing angle between the two objects [rad]
% V1 - orbital velocity of the chaser
% e, i, w, om - eccentricity, incination, RAAN, argument of periapsis

%% Outputs

% DeltaV1 - first impulse of the Hohmann transfer
% DeltaV2 - second impulse of the Hohmann transfer
% ToF - time of flight of the transfer
% phi - phase angle at which the transfer starts
% t_phasing - time to wait before starting the transfer

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 13/05/2025
% Update date:

addpath("Functions\General_Functions")
addpath("Functions\Function_Part1")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hohmann Transfert - Analitical 

% Computing the Hohmann transfer deltaV using the analitical formulation
[DeltaV1,DeltaV2, DeltaVtot, a_h, e_h, ToF] = HohmannTransferAnalitical(r1,r2, mu_Earth);


fprintf('DeltaV1: %.4f km/s\n', DeltaV1);
fprintf('DeltaV2: %.4f km/s\n', DeltaV2);
fprintf('Total DeltaV: %.4f km/s\n', DeltaVtot);
fprintf('Time of Flight: %.1f s (%.2f minutes)\n', ToF, ToF/60);

%% Phasing

% Mean motion (= angular velocity) for the two bodies
n2=sqrt(mu_Earth/r2^3);
n1=sqrt(mu_Earth/r1^3);
%%%%%%%
%t_phasing=(n2*ToF-pi-phi0_rad)/(n1-n2);
phi = sqrt(mu_Earth/r2^3)*pi*sqrt((r1+r2)^3/8/mu_Earth) - pi;
%phi=phi0_rad-(n2-n1)*t_phasing;
t_phasing = (phi0_rad-phi)/(n2-n1);

fprintf('Phasing angle: %.4f rad\n', phi);
fprintf('Waiting time: %.1f s (%.2f minutes)\n \n', t_phasing, t_phasing/60);

%% Plot of the analyitc orbits

[Htransfer, orbit1, orbit2] = plot_Analytical_orbit(a_h,e_h, r1, r2, n1, n2,t_phasing,phi);


%% Hohmann Transfer - Numerical 

% Assumption: Hohmann transfer starts at pi cartesian
f_start = pi;

% Angular momentum at the starting point of the Homann transfer
h_start_H = (V1+DeltaV1)*r1;

% Numerical solution of the Hohmann transfer
[t,r,r_dot,f] = HohmannTransferNumerical(r1,e,i, om, w, f_start, mu_Earth, ToF,h_start_H, DeltaV1);




%% Plot numerical obrit

plot_numerical(Htransfer,r, a_h,e_h,f)

%% Constant of motion

[error_h_rel,error_e_rel,error_E_rel] = constant_of_motion(r,r_dot, t,r1,h_start_H,mu_Earth, e_h, V1 , DeltaV1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hohmann transfer - Numerical v2
f_start = 0;
[t,r,r_dot,~] = HohmannTransferNumerical_v2(r1,e,i, om, w, f_start, mu_Earth, ToF,h_start_H, DeltaV1,DeltaV2, t_phasing);

figure
axis equal; hold on; grid on;

% Plot complete circular orbits
plot(orbit1(1,:), orbit1(2,:), 'b', 'LineWidth', 0.5);
plot(orbit2(1,:), orbit2(2,:), 'r', 'LineWidth', 0.5);
plot(0, 0, 'bo', 'MarkerSize', 20, 'MarkerFaceColor','k','HandleVisibility','off');

% Plot the numerical Hohman transfer
plot(r(1,:),r(2,:), 'k-', 'LineWidth', 2);
plot(r(1,1),r(2,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b','HandleVisibility','off');
plot(r(1,end),r(2,end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r','HandleVisibility','off');
title('Numerical and analitical Hohmann transfer')
legend("Starting orbit","Final Orbit", "Numerical tranfer");
xlabel('X [km]'); ylabel('Y [km]');

%% Constant of motion :

% 1. Angular momentum h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_vec = cross(r, r_dot);
h_value = vecnorm(h_vec);

% 2. Eccentricity vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
e_vec = cross(r_dot, h_vec)/mu_Earth - r./vecnorm(r);
e_value = vecnorm(e_vec);

% 3. Total energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total orbital energy of the Hohmann transfer (numerical)
E = 0.5 * vecnorm(r_dot).^2 - mu_Earth ./ vecnorm(r);

%% Plot of the results 

figure;
subplot(3,1,1);
plot(t, h_value, 'r', 'LineWidth', 0.5); hold on; grid on;
title('Angular momentum'); ylabel('h $ [\mathrm{km^2}/\mathrm{s}]$', 'Interpreter', 'latex');
xlabel('Time [s]');

subplot(3,1,2);
plot(t, e_value, 'g', 'LineWidth', 0.5); hold on; grid on;
title('Eccentricity vector'); ylabel('$e [-]$', 'Interpreter', 'latex');
xlabel('Time [s]');

subplot(3,1,3);
plot(t, E, 'b', 'LineWidth', 0.5); hold on; grid on;
title('Total orbital energy'); ylabel('$\epsilon$', 'Interpreter', 'latex');
xlabel('Time [s]');


end