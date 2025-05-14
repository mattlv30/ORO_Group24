clear all; clc; 

% r = a + (b-a).*rand(n,1)
rH_err = [500 + (1000-500).*rand(2,1);0]./1000; % errors between 500 and 1000 in the plane [m]-->[km]
rH_dot_err = [0.1 + (1-0.1).*rand(2,1);0]./1000; % errors between 0.1 and 1in the plane [m/s]-->[km/s]