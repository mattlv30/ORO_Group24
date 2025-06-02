%% Solve KR2BP (specific for part 3)

%% Input 
% mu, h =sqrt(mu*a*(1-e^2), intial condition r0 r_dot0

%% Output time, r, r_dot
% informations KR2BP [r_dot_dot]= -mu*[r]/r^3 , [] a vector , r norm

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 30/05/2025
% Update date:

function [t,r,r_dot] = solve2bodyK_part3(mu,cond_i,tspan)
% z=[x y z vx vy vz f]
% standard substitution
ff=@(t,z) [z(4);z(5);z(6);
    -mu*z(1)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
    -mu*z(2)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
    -mu*z(3)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3];
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','off'); % to change tol

[t, z] = ode113(ff,tspan,cond_i,opts);
r=z(:,1:3)'; % every column has the cartesian coordinates of r at that time
r_dot=z(:,4:6)';  % every column has the cart coord of r_dot at that time
end