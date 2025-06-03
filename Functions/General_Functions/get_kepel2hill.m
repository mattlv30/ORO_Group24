%% Orbital element to Hill reference frame  
% Description from Kepler orbital elements to state vector Hill frame

%% Input 
% semimajor axis a, eccentricity e, true anomaly f, grav parameter mu

%% Output 
% state vector in Hill frame rH rH_dot (derivative) [rH, rH_dot]

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 09/05/2025
% Update date:

function [rH,rH_dot] = get_kepel2hill(a,e,f,mu)
rH= a*(1-e^2)/(1+e*cos(f)) *[1; 0; 0];
rH_dot= sqrt(mu/a /(1-e^2))*[e*sin(f);1+e*cos(f);0];
end