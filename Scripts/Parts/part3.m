%% Part 3
% This function solves all the request of the third part of the project

function [] = part3(rho,CD,A,m,ToF_min,r2,mu_Earth,rH_err,rH_dot_err,deltaV1,deltaV2)

%% Inputs
% rho - desnity
% CD - drag coefficient
% A- cross sectional area
% m - satellite mass
% ToF_min - ToF of the minimum Delta V trajcetory
% r2 - radius of the target orbit
% mu_Earth - gravitational pararmeter of the primary (Earth)
% rH_err - random error in position
% rH_dot_err - random error in velocity
% deltaV1 - first DeltaV
% deltaV2 - second DeltaV

%% Outputs

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 13/05/2025
% Update date:

addpath("Functions\General_Functions")
addpath("Functions\Function_Part2")
addpath("Functions\Function_Part3")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n2=sqrt(mu_Earth/r2^3);
% fixed tspan to confront different propagations
tspan=linspace(0, ToF_min,500); % adimensional time

% get initial condition target in Equatorial
[rH_target,rH_dot_target] = get_kepel2hill(r2,0,0,mu_Earth);
[rP_target,rP_dot_target] = get_hill2perifocal(rH_target,rH_dot_target,0);
[rE_target,rE_dot_target] = get_perifocal2equatorial(rP_target,rP_dot_target,0,0,0);

%% check rendezvous in Equatorial restricted 2bodyK without perturbation

% adimensional propagation
[~,r_t,~] = solve2bodyK_part3(1,[rE_target./r2;rE_dot_target./r2./n2],tspan);
r_t=r_t*r2; % dimensional position target

% initial condition chaser in Equatorial
r_c0=rH_err+[r2;0;0]; 
r_dot_c0=rH_dot_err+[0;n2*r2;0]+deltaV1*r2*n2+cross([0;0;n2],rH_err);

% adimensional propagation
[t,r_c,~] = solve2bodyK_part3(1,[r_c0/r2;r_dot_c0/r2/n2],tspan);
r_c=r_c*r2; % dimensional position chaser
t=t/n2; % dimensional time = tt

figure;
hold on
plot(r_t(1,:),r_t(2,:),"r")
plot([r_t(1,1),r_c(1,1)],[r_t(2,1),r_c(2,1)],"k^","LineWidth",1)
plot([r_t(1,end),r_c(1,end)],[r_t(2,end),r_c(2,end)],"g*","LineWidth",1)
plot(r_c(1,:),r_c(2,:),"b")
axis equal
title("Rendezvous without drag")
xlabel("x [km]")
ylabel("y [km]")
legend("Target","Starting points","Ending points","Chaser")
hold off

figure;
plot(t,vecnorm(r_t-r_c,2,1),"r")
title("Distance target - chaser without drag")
xlabel("t [s]")
ylabel("[km]")

disp(" ==================== Part 3 - Perturbations ==================== ")
disp(" ")

fprintf('Final position offset target - chaser without drag (equatorial CCS): %.4f m \n', norm(r_t(:,end)-r_c(:,end),2)*1000);

%% rendezvous in Equatorial with drag on chaser
% drag adimensionalizzata
% rho assumed constant since h~350km during the maneuver
coeff_ad=1/2*rho*CD*A/m *r2;
% z = [x;y;z;vx;vy;vz]
% standard substitution
d_dt_ad=@(t,z) [
    z(4);
    z(5);
    z(6);
    % adimensionalization mu-->1
    -1*z(1)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3-coeff_ad*sqrt(z(4)^2+z(5)^2+z(6)^2)*z(4);
    -1*z(2)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3-coeff_ad*sqrt(z(4)^2+z(5)^2+z(6)^2)*z(5);
    -1*z(3)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3-coeff_ad*sqrt(z(4)^2+z(5)^2+z(6)^2)*z(6)];

opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','off');
[tc_drag,sol_chaser_drag]=ode113(d_dt_ad,tspan,[r_c0/r2;r_dot_c0/r2/n2],opts);
r_c_drag=sol_chaser_drag(:,1:3)'*r2; % dimensional chaser position vector
% r_c_dot_drag=sol_chaser_drag(:,4:6)'*r2; % dimensional chaser velocity vector
tc_drag=tc_drag/n2; % dimensional time = tt from target numerical integration



figure;
plot(tc_drag,vecnorm(r_c_drag-r_c,2,1),"b") % position difference chaser perturbated/non
title("Distance chaser with drag from its unperturbated orbit")
xlabel("t [s]")
ylabel("[km]")

figure;
hold on
plot(r_t(1,:),r_t(2,:),"r")
plot([r_t(1,1),r_c_drag(1,1)],[r_t(2,1),r_c_drag(2,1)],"k^","LineWidth",1)
plot([r_t(1,end),r_c_drag(1,end)],[r_t(2,end),r_c_drag(2,end)],"g*","LineWidth",1)
plot(r_c_drag(1,:),r_c_drag(2,:),"b")
axis equal
title("Rendezvous with drag on chaser")
xlabel("x [km]")
ylabel("y [km]")
legend("Target","Starting points","Ending points","Chaser")
hold off

figure;
plot(t,vecnorm(r_t-r_c_drag,2,1),"r")
title("Distance target - chaser with drag")
xlabel("t [s]")
ylabel("[km]")

fprintf('Final position offset target - chaser with drag (equatorial CCS): %.4f m \n', norm(r_t(:,end)-r_c_drag(:,end),2)*1000);

%% rendezvous with drag in the relative CCS with HCW assumptions
[~,csi,eta,~,v_csi,v_eta,~] = solveHCW_drag(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV1],tspan);

figure;
hold on
axis equal
plot(csi,eta,"b")
plot(csi(end),eta(end),"r*")
title("Relative position chaser with drag")
xlabel("\xi [-]")
ylabel("\eta [-]")
hold off

figure;
hold on
axis equal
plot(v_csi,v_eta,"b",'HandleVisibility','off')
plot(v_csi(end)+deltaV2(1),v_eta(end)+deltaV2(2),"r*")
title("Relative velocity chaser with drag")
xlabel("v_{\xi} [-]")
ylabel("v_{\eta} [-]")
legend("V after 2° impulse")
hold off

fprintf('Final position offset target - chaser with drag relative CCS: %.4f m\n', norm([csi(end);eta(end)],2)*r2*1000);
fprintf('Final velocity offset target - chaser with drag relative CCS: %.4f m/s \n', norm([v_csi(end)+deltaV2(1);v_eta(end)+deltaV2(2)],2)*r2*1000);

disp(" ")

end

% fprintf('Final velocity offset target - chaser with drag relative CCS: %.4f [m/s]\n', ([v_csi(end);v_eta(end)]+deltaV2([1 2]))*r2*1000);

% %% PLOT IDEAL AND DRAG TRAJECTORY TOGETHER
% [~,csi_id,~,eta_id,~,~,~] = solveHCW([rH_err(1)/r2;rH_dot_err(1)/r2/n2+deltaV1(1);rH_err(2)/r2;rH_dot_err(2)/r2/n2+deltaV1(2);rH_err(3)/r2;rH_dot_err(3)/r2/n2+deltaV1(3)],tspan);
% figure;
% hold on
% axis equal
% plot(csi_id,eta_id,"b")
% plot(csi_id(end),eta_id(end),"b*")
% plot(csi,eta,"r")
% plot(csi(end),eta(end),"r*")
% title("Relative position chaser")
% xlabel("\xi [-]")
% ylabel("\eta [-]")
% hold off

