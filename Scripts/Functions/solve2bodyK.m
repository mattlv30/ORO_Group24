%% HohmannTransferNumerical
% This function numerically solves the Kepler restricted 2 body problem. It
% also computes the numerical evolution of the true anomaly of the
% satellite

%% Inputs
% mu - gravitational constant of the primary 
% h - angular momentum computed at the starting point
% cond_i - initial condition: 3 positions, 3 velocities, 1 angular momentum
% tspan - time span (e.g. the orbital oeriod)

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
function [t,r,r_dot,f] = solve2bodyK(mu,h,cond_i,tspan)

% Standard substitution of a system of ODE: z=[x y z vx vy vz f]

% Function used (legend: [r] a vector , r norm) :  
% KR2BP: [r_dot_dot] = -mu*[r]/r^3  
% True anomaly: f_dot = h/r^2

ff=@(t,z) [
    z(4);
    z(5);
    z(6);
   -mu*z(1)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
   -mu*z(2)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
   -mu*z(3)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
    h/(z(1)^2+z(2)^2+z(3)^2) % in the KR2BP h=cost
    ]; 

% Change in the tolerance to improve numerical integration
opts = odeset('Reltol',1e-8,'AbsTol',1e-9,'Stats','off'); 

% Numerical solution using ode113
[t, z] = ode113(ff,tspan,cond_i, opts);
r=z(:,1:3)'; % every column has the cartesian coordinates of r at that time
r_dot=z(:,4:6)';  % every column has the cart coord of r_dot at that time
f=z(:,7)'; % the column vector containing the true anomaly

end