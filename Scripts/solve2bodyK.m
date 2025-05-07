% Description 
% input mu, h =sqrt(mu*a*(1-e^2), intial condition r0 r_dot0
% output time, r, r_dot
% informations KR2BP [r_dot_dot]= -mu*[r]/r^3 , [] a vector , r norm
% author CI
% creation date 18/03/2025
% update date

function [t,r,r_dot,f] = solve2bodyK(mu,h,cond_i,tspan)
% z=[x y z vx vy vz f]
% standard substitution
ff=@(t,z) [z(4);z(5);z(6);
    -mu*z(1)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
    -mu*z(2)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
    -mu*z(3)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
    h/(z(1)^2+z(2)^2+z(3)^2)]; 
% opts = odeset('Reltol',1e-8,'AbsTol',1e-9,'Stats','on'); % to change tol

[t, z] = ode113(ff,tspan,cond_i);
r=z(:,1:3)'; % every column has the cartesian coordinates of r at that time
r_dot=z(:,4:6)';  % every column has the cart coord of r_dot at that time
f=z(:,7)';
end