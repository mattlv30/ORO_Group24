clc; clear; close all;

%% Constants
mu = 398600.4418;          % [km^3/s^2] Earth’s gravitational parameter
Re = 6378.1;               % [km] Earth radius
h1 = 600;                  % [km] initial orbit altitude (S/C)
h2 = 300;                  % [km] target orbit altitude (Debris)
R1 = Re + h1;
R2 = Re + h2;

%% Velocities and Transfer Orbit
v1 = sqrt(mu / R1);        % Velocity on initial orbit
v2 = sqrt(mu / R2);        % Velocity on target orbit

a_trans = (R1 + R2)/2;                         % Semi-major axis of transfer ellipse
v_trans1 = sqrt(2*mu*R2 / (R1*(R1+R2)));       % Velocity after first impulse
v_trans2 = sqrt(2*mu*R1 / (R2*(R1+R2)));       % Velocity before second impulse

deltaV1 = v_trans1 - v1;
deltaV2 = v2 - v_trans2;
deltaV_total = deltaV1 + deltaV2;

fprintf('DeltaV1: %.4f km/s\n', deltaV1);
fprintf('DeltaV2: %.4f km/s\n', deltaV2);
fprintf('Total DeltaV: %.4f km/s\n', deltaV_total);

%% Time of Flight
ToF = pi * sqrt(a_trans^3 / mu);  % Half orbit period
fprintf('Time of Flight: %.1f s (%.2f minutes)\n', ToF, ToF/60);

%% Plotting the Orbits
theta = linspace(0, 2*pi, 500);
circ1 = R1 * [cos(theta); sin(theta)];
circ2 = R2 * [cos(theta); sin(theta)];

% Transfer ellipse
theta_trans = linspace(0, pi, 250);
r_trans = (R1 * R2) ./ ((R2 - R1)/2 * cos(theta_trans) + (R1 + R2)/2);
x_trans = r_trans .* cos(theta_trans);
y_trans = r_trans .* sin(theta_trans);

figure;
hold on; grid on; axis equal;
plot(circ1(1,:), circ1(2,:), 'b--', 'LineWidth', 1.5);
plot(circ2(1,:), circ2(2,:), 'r--', 'LineWidth', 1.5);
plot(x_trans, y_trans, 'k-', 'LineWidth', 2);
plot(R1, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b');
plot(R2*cos(pi), R2*sin(pi), 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
legend('Initial Orbit', 'Target Orbit', 'Hohmann Transfer', 'Start Point', 'End Point');
title('Hohmann Transfer Maneuver');
xlabel('X [km]'); ylabel('Y [km]');

%% Numerical Propagation to Check Constants of Motion
twoBody = @(t, y) [y(4:6); -mu * y(1:3) / norm(y(1:3))^3];
r0 = [R1; 0; 0];
v0 = [0; v_trans1; 0];
y0 = [r0; v0];

[t_out, y_out] = ode45(twoBody, [0 ToF], y0);

r_mag = vecnorm(y_out(:,1:3),2,2);
v_mag = vecnorm(y_out(:,4:6),2,2);
E = 0.5 * v_mag.^2 - mu ./ r_mag;
h_vec = cross(y_out(:,1:3), y_out(:,4:6), 2);
h_mag = vecnorm(h_vec,2,2);

figure;
subplot(2,1,1);
plot(t_out, h_mag); grid on;
title('Angular Momentum Magnitude'); ylabel('|h| [km^2/s]');

subplot(2,1,2);
plot(t_out, E); grid on;
title('Specific Orbital Energy'); ylabel('Energy [km^2/s^2]');
xlabel('Time [s]');