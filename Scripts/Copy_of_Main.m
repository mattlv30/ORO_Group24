%% ORO Project Group 24

clear; clc; close all; 

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

% Orbital Velocity 
V1 = sqrt(mu_Earth/r1); %[km/s]
V2 = sqrt(mu_Earth/r2); %[km/s]
n2=sqrt(mu_Earth/r2^3);

%% part 2
[rH_err,rH_dot_err,Delta_T,deltaV1,deltaV2] = part2(r2,n2);
%% part3 be careful unit of measurement
rho=10^(-12)*10^9; % atmospheric density [kg/m^3] ---> [kg/km^3]
CD=2.2; % drag coefficient
A=2/10^6; % cross-sectional area sat [m^2] ---> [km^2]
m=325; % sat mass [kg]

tspan=linspace(0, Delta_T,500); %tspan=[0 Delta_T];
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');

[rH_target,rH_dot_target] = get_kepel2hill(r2,0,0,mu_Earth);
[rP_target,rP_dot_target] = get_hill2perifocal(rH_target,rH_dot_target,0);
[rE_target,rE_dot_target] = get_perifocal2equatorial(rP_target,rP_dot_target,0,0,0);

cond_i_target=[rE_target(1)/r2;rE_dot_target(1)/r2/n2;rE_target(2)/r2;rE_dot_target(2)/r2/n2;rE_target(3)/r2;rE_dot_target(3)/r2/n2];
[tt,r_t,r_dot_t] = Copy_of_solve2bodyK(1,[rE_target./r2;rE_dot_target./r2./n2],tspan);
r_t=r_t*r2;
tt=tt/n2;

r_c0=rH_err+[r2;0;0]; r_c0=r_c0/r2;
r_dot_c0=rH_dot_err+[0;n2*r2;0]+deltaV1*r2*n2+cross([0;0;n2],rH_err); r_dot_c0=r_dot_c0/r2/n2;
[t,r,r_dot] = Copy_of_solve2bodyK(1,[r_c0;r_dot_c0],tspan);
r=r*r2; % posizione ECI chaser
t=t/n2;

% rendezvous con target cerchio analitico
figure;
hold on
% plot(r_t(1,:),r_t(2,:),"r")
% plot(r_t(1,end),r_t(2,end),"ro","LineWidth",3)
plot(r2*cos(linspace(0,2*pi,1000)),r2*sin(linspace(0,2*pi,1000)),"r")
plot(r(1,:),r(2,:),"b")
plot(r(1,1),r(2,1),"bo")
plot(r(1,end),r(2,end),"b*")
axis equal
hold off

% rendezvous in equatoriale con target da integrazione
figure;
hold on
plot(r_t(1,:),r_t(2,:),"r")
plot(r_t(1,1),r_t(2,1),"ro","LineWidth",3)
plot(r_t(1,end),r_t(2,end),"r*","LineWidth",3)
plot(r(1,:),r(2,:),"b")
plot(r(1,1),r(2,1),"bo")
plot(r(1,end),r(2,end),"b*")
axis equal
hold off

%drag adimensionalizzata
% z = [x;vx;y;vy;z;vz]
coeff_ad=1/2*rho*CD*A/m *r2;
d_dt_ad=@(t,z) [
    z(2);
    -1*z(1)/sqrt(z(1)^2+z(3)^2+z(5)^2)^3-coeff_ad*sqrt(z(2)^2+z(4)^2+z(6)^2)*z(2);
    z(4);
    -1*z(3)/sqrt(z(1)^2+z(3)^2+z(5)^2)^3-coeff_ad*sqrt(z(2)^2+z(4)^2+z(6)^2)*z(4);
    z(6);
    -1*z(5)/sqrt(z(1)^2+z(3)^2+z(5)^2)^3-coeff_ad*sqrt(z(2)^2+z(4)^2+z(6)^2)*z(6)];
% dr=0

% soluzioni kepler perturbate adimensionalizzate
[tt_drag,sol_target]=ode113(d_dt_ad,tspan,cond_i_target,opts);
[tc_drag,sol_chaser]=ode113(d_dt_ad,tspan,[r_c0(1);r_dot_c0(1);r_c0(2);r_dot_c0(2);r_c0(3);r_dot_c0(3)],opts);
figure;
axis equal
plot(sol_target(:,1)*r2,sol_target(:,3)*r2,"r")
hold on
plot(sol_chaser(:,1)*r2,sol_chaser(:,3)*r2,"b")
hold off

figure; % length needs to be equal
subplot(2,1,1);
plot(t,vecnorm(r_t-r,2),"r") % distanza target chaser non perturbato
subplot(2,1,2);
plot(tt_drag/n2,vecnorm(sol_chaser(:,[1 3])-sol_target(:,[1 3]),2,2)*r2,"b") % distanza target chaser perturbato

figure; % length needs to be equal
subplot(2,1,1);
plot(t,vecnorm(r_t([1 2],:)'-sol_target(:,[1 3])*r2,2,2),"r") % differenza target non perturbata e con drag
subplot(2,1,2);
plot(tt_drag/n2,vecnorm(sol_chaser(:,[1 3])*r2-r([1 2],:)',2,2),"b") % differenza chaser non e perturbato

r_t([1 2],end)
r([1 2],end)
r([1 2],end)-r_t([1 2],end)

sol_chaser(end,[1 3])*r2
sol_target(end,[1 3])*r2
sol_chaser(end,[1 3])*r2-sol_target(end,[1 3])*r2

% chiedere se drag da applicare solo al chaser
% caso con drag solo al chaser

figure;
axis equal
plot(r_t(1,:),r_t(2,:),"r")
hold on
plot(sol_chaser(:,1)*r2,sol_chaser(:,3)*r2,"b")
hold off

figure; % length needs to be equal
plot(t,vecnorm(r_t([1 2],:)'-sol_chaser(:,[1 3])*r2,2,2),"r") % rendezvou chaser perturbato, target non
sol_chaser(end,[1 3])*r2-r_t([1 2],end)'

%% implementare hill con chaser v fissa