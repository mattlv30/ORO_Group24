%% Part 2
% This function solves all the request of the second part of the project

function [rH_err,rH_dot_err,Delta_T,deltaV1,deltaV2] = part2(r2,n2)

%% Inputs
% r2 - radius of the target orbit
% n2 - mean motion of the target object

%% Outputs

% rH_err - random error in position
% rH_dot_err - random error in velocity
% Delta_T - Delta T (ToF) for min deltaV
% deltaV1 - first DeltaV
% deltaV2 - second DeltaV

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 13/05/2025
% Update date:

addpath("Functions\")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Parts\'Random error_part2'\rH_err.mat rH_err
load Parts\'Random error_part2'\rH_dot_err.mat rH_dot_err


% compute for different Delta T time of flight (f_)
% to adimentionalize it is used L=r2 n=n2 of the target
% using L=R for normalization distance
rH_err_ad=rH_err./r2;
% 1 /L /n for normalization velocity
rH_dot_err_ad=rH_dot_err./r2./n2;

ff0=0; % initial epoch
ffend=2*2*pi; % two orbital periods of the target normalized
N=1000; % number of evaluations
ff=linspace(ff0,ffend,N);
ff=ff(2:end); % eliminate 0 from the possible ToF
f0=0; % initial time of the maneuver

% adimetionalized Hill components
csi0=rH_err_ad(1);
eta0=rH_err_ad(2);
zeta0=rH_err_ad(3);
v_csi0=rH_dot_err_ad(1);
v_eta0=rH_dot_err_ad(2);
v_zeta0=rH_dot_err_ad(3);

deltaV1=zeros(N-1,3);
deltaV2=zeros(N-1,3);
deltaVtot=zeros(N-1,1);

for ii=1:length(ff)
    f_=ff(ii); % priori chosen maneuver duration time in normalized unit

    % first impulse to get in the target position [0 0 0] Hill RF
    % vplus velocity after impulse
    vplus_csi0=(-2*eta0+(2*eta0-3*f_*csi0)*cos(f_)+4*csi0*sin(f_))/(3*f_*sin(f_)+8*cos(f_)-8);
    vplus_eta0=(14*csi0-14*csi0*cos(f_)+(eta0-6*f_*csi0)*sin(f_))/(3*f_*sin(f_)+8*cos(f_)-8);
    vplus_zeta0=-cos(f_)*zeta0/sin(f_);
    delta1=[vplus_csi0-v_csi0; vplus_eta0 - v_eta0; vplus_zeta0-v_zeta0];
    
    [phi_p, phi_c] = analytical_solHCW(f0); % matrix
    Ap=phi_p(f_);
    sol_csi_eta=Ap*[csi0;vplus_csi0;eta0;vplus_eta0];

    % the position components are commented since they are not used
    % csi=sol_csi_eta(1,:);
    % eta=sol_csi_eta(3,:);
    v_csi=sol_csi_eta(2,:);
    v_eta=sol_csi_eta(4,:);
    Ac=phi_c(f_);

    % the zeta components wouldn't be needed since they are 0
    sol_zeta=Ac*[zeta0;vplus_zeta0];
    % zeta=sol_zeta(1,:);
    v_zeta=sol_zeta(2,:);
    
    % second impulse to match the target's velocity
    delta2=-[v_csi(end),v_eta(end),v_zeta(end)]';
    
    % total deltaV required for rendezvous single directional thruster is used
    deltaV1(ii,:)=delta1;
    deltaV2(ii,:)=delta2;
    deltaVtot(ii)= norm(delta1,2)+norm(delta2,2);

end

[min_dV, idx]=min(deltaVtot); % find minimum total impulse needed
Delta_T=ff(idx); % select Delta T (ToF) for min deltaV
deltaV1=deltaV1(idx,:)';
deltaV2=deltaV2(idx,:)';
fprintf('Minimum total impulse: %.4f [km/s]\n', min_dV*r2);
fprintf('ToF for minimum DeltaV: %.4f [s] (%.2f minutes)\n', Delta_T/n2, Delta_T/n2/60);

figure;
plot(ff,deltaVtot,Color="b")
title("Values of total impulse for different ToF")
ylim([0 0.2])
xlabel("ToF [-]")
ylabel("\Delta V [-]")

% numerical solution for chosen impulse and ToF
vplus_csi0=v_csi0+deltaV1(1);
vplus_eta0=v_eta0+deltaV1(2);
vplus_zeta0=v_zeta0+deltaV1(3);
[t,csi,v_csi,eta,v_eta,~,~]=solveHCW([csi0;vplus_csi0;eta0;vplus_eta0;zeta0;vplus_zeta0],[0 Delta_T]);

% difference between analytical and numerical solutions
err_csi=zeros(length(t),1);
err_v_csi=zeros(length(t),1);
err_eta=zeros(length(t),1);
err_v_eta=zeros(length(t),1);
for kk=1:length(t)
    Ap=phi_p(t(kk));
    sol_csi_eta=Ap*[csi0;vplus_csi0;eta0;vplus_eta0];
    csi_an=sol_csi_eta(1,:);
    eta_an=sol_csi_eta(3,:);
    v_csi_an=sol_csi_eta(2,:);
    v_eta_an=sol_csi_eta(4,:);

    % zeta components are 0
    % Ac=phi_c(t(kk));
    % sol_zeta=Ac*[zeta0;vplus_zeta0];
    % zeta_an=sol_zeta(1,:);
    % v_zeta_an=sol_zeta(2,:);

    err_csi(kk)=abs(csi(kk)-csi_an);
    err_v_csi(kk)=abs(v_csi(kk)-v_csi_an);
    err_eta(kk)=abs(eta(kk)-eta_an);
    err_v_eta(kk)=abs(v_eta(kk)-v_eta_an);
end

% plot dimensional errors
figure;
plot(t*n2,err_csi*r2)
title("Error postion h_1")
xlabel("t [s]")
ylabel("[km]")

figure;
plot(t*n2,err_v_csi*r2*n2)
title("Error position h_2")
xlabel("t [s]")
ylabel("[km]")

figure;
plot(t*n2,err_eta*r2)
title("Error velocity h_1")
xlabel("t [s]")
ylabel("[km/s]")

figure;
plot(t*n2,err_v_eta*r2*n2)
title("Error velocity h_2")
xlabel("t [s]")
ylabel("[km/s]")

% rendezvous trajectory in Hill CCS
figure;
plot(csi,eta,"b")
title("HWC solution position")
xlabel("\xi [-]")
ylabel("\eta [-]")

% velocity evolution in Hill CCS
figure;
hold on
plot(v_csi,v_eta,'HandleVisibility','off')
plot(v_csi(end)+deltaV2(1),v_eta(end)+deltaV2(2),"r*")
hold off
title("HWC solution velocity")
xlabel("v_{\xi} [-]")
ylabel("v_{\eta} [-]")
legend("V after 2° impulse")

end