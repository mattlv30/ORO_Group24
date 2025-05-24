%% HohmannTransferNumerical
% This function numerically computes the Hohmann transfer.

% It perform 3 changes of coordinates going from the orbital elements to
% the perifocal reference frame. Then the impulsive
% velocity needed for the Hohmann transfer is added to the orbital
% velocity. Finally, the numerical solution of the problem is computed
% using a numerical solver

%% Inputs
% a, e, i, om, w, f - classical orbital elements (semi major axis,
% eccentricity, inclination, RAAN, argument of periapsis, ture anomaly) at
% the starting point of the transfer
% mu - gravitational constant of the primary 
% ToF - time of flight of the transfer
% h - angular momentum computed at the starting point of the transfer
% DeltaV1 - deltaV to start the manouvre computed with the Hohmann trnasfer
% analitical formula

%% Outputs

% t - time evaluation point from the numerical solver (ode113)
% r - vector position in the ECI reference frame [x; y; z]
% r_dot - velocity vector in the ECI reference frame [Vx; Vy; Vz]
% f - ture anomaly

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 09/05/2025
% Update date:

%% Function 
function [t,r,r_dot,f] = HohmannTransferNumerical(a,e,i, om, w, f, mu, ToF,h, DeltaV1)

% To cartesian Hill reference frame
[rH0, rH_dot0] = get_kepel2hill(a, e, f, mu); 

% The impulseive manouvre for the hohmann transfer is added here as a deltaV
% The impules is assumed to be tangent to the orbit (theta direction in the
% Hill reference frame)

[rP0, rP_dot0] = get_hill2perifocal(rH0,rH_dot0 + [0 DeltaV1 0]',f); % to perifocal
[rE0, rE_dot0] = get_perifocal2equatorial(rP0,rP_dot0,om,i,w); % to equatorial

% Definition of the time span (time of flight of the Homann transfer)
tspan = linspace(0, ToF)';

% Definition of the 7 initial condition: 3 position, 3 velocities, 1
% angular momentum 
cond_i = [rE0; rE_dot0; f];

% Integration of the problem using a function that solves the Kepler
% restricted 2 body problem with a numerical integration (ode113)
[t,r,r_dot,f]=solve2bodyK(mu,h,cond_i,tspan);


end