% part 2
% decidere se tenere variabili non usate
% migliorare i plot

function [rH_err,rH_dot_err,Delta_T,deltaV1,deltaV2] = part2(r2,n2)
% r = a + (b-a).*rand(n,1)
% randn normal distribution does not have min or max values

% check unit of measurement
rH_err = [500 + (1000-500).*rand(2,1);0]./1000; % errors between 500 and 1000 in the plane [m]-->[km]
rH_dot_err = [0.1 + (1-0.1).*rand(2,1);0]./1000; % [m/s]--->[km/s]
% compute for different Delta T time of flight (f_)
% for normalization it is used L=r2 n=n2 target
% using L=R for normalization distance
rH_err_ad=rH_err./r2;
% 1 /L /n for normalization velocity
rH_dot_err_ad=rH_dot_err./r2./n2;
ff0=0;
ffend=2*2*pi; % two orbital periods of the target normalized
N=1000;
ff=linspace(ff0,ffend,N);
f0=0; % initial time of the maneuver

csi0=rH_err_ad(1);
eta0=rH_err_ad(2);
zeta0=rH_err_ad(3);
v_csi0=rH_dot_err_ad(1);
v_eta0=rH_dot_err_ad(2);
v_zeta0=rH_dot_err_ad(3);

deltaV1=zeros(N,3);
deltaV2=zeros(N,3);
deltaVtot=zeros(N,1);

for ii=1:length(ff)
    f_=ff(ii); % priori chosen maneuver duration time in normalized unit
    % first impulse to get in the target position [0 0 0] Hill RF
    % vplus velocity after impulse
    vplus_csi0=(-2*eta0+(2*eta0-3*f_*csi0)*cos(f_)+4*csi0*sin(f_))/(3*f_*sin(f_)+8*cos(f_)-8);
    vplus_eta0=(14*csi0-14*csi0*cos(f_)+(eta0-6*f_*csi0)*sin(f_))/(3*f_*sin(f_)+8*cos(f_)-8);
    vplus_zeta0=-cos(f_)*zeta0/sin(f_);
    delta1=[vplus_csi0-v_csi0; vplus_eta0 - v_eta0; vplus_zeta0-v_zeta0];
    
    [phi_p, phi_c] = analytical_solHCW(f0);
    Ap=phi_p(f_);
    sol_csi_eta=Ap*[csi0;vplus_csi0;eta0;vplus_eta0];
    % csi=sol_csi_eta(1,:);
    % eta=sol_csi_eta(3,:);
    v_csi=sol_csi_eta(2,:);
    v_eta=sol_csi_eta(4,:);
    Ac=phi_c(f_);
    sol_zeta=Ac*[zeta0;vplus_zeta0];
    % zeta=sol_zeta(1,:);
    v_zeta=sol_zeta(2,:);
    
    % vplus to final velocity befor second impulse <-- gravitational effects
    % second impulse
    delta2=-[v_csi(end),v_eta(end),v_zeta(end)]';
    
    % total deltaV required for rendezvous
    deltaV1(ii,:)=delta1;
    deltaV2(ii,:)=delta2;
    deltaVtot(ii)= norm(delta1,2)+norm(delta2,2);

end

[min_dV, idx]=min(deltaVtot);
Delta_T=ff(idx); % select Delta T for min deltaV
deltaV1=deltaV1(idx,:)';
deltaV2=deltaV2(idx,:)';
fprintf('min_dV: %.4f km/s\n', min_dV);
fprintf('ToF for min_dV: %.4f km/s\n', Delta_T);

figure;
plot(ff,deltaVtot,Color="b")
ylim([0 0.5])
% how to treat first NaN (impossibile 0 seconds maneuver)

vplus_csi0=v_csi0+deltaV1(1);
vplus_eta0=v_eta0+deltaV1(2);
vplus_zeta0=v_zeta0+deltaV1(3);
[t,csi,v_csi,eta,v_eta,~,~]=solveHCW([csi0;vplus_csi0;eta0;vplus_eta0;zeta0;vplus_zeta0],[0 Delta_T]);

err_csi=zeros(length(t),1);
err_v_csi=zeros(length(t),1);
err_eta=zeros(length(t),1);
err_v_eta=zeros(length(t),1);
for kk=1:length(t)
    [phi_p, ~] = analytical_solHCW(f0);
    Ap=phi_p(t(kk));
    sol_csi_eta=Ap*[csi0;vplus_csi0;eta0;vplus_eta0];
    csi_an=sol_csi_eta(1,:);
    eta_an=sol_csi_eta(3,:);
    v_csi_an=sol_csi_eta(2,:);
    v_eta_an=sol_csi_eta(4,:);
    % Ac=phi_c(t(kk));
    % sol_zeta=Ac*[zeta0;vplus_zeta0];
    % zeta_an=sol_zeta(1,:);
    % v_zeta_an=sol_zeta(2,:);

    err_csi(kk)=abs(csi(kk)-csi_an);
    err_v_csi(kk)=abs(v_csi(kk)-v_csi_an)/abs(v_csi_an);
    err_eta(kk)=abs(eta(kk)-eta_an);
    err_v_eta(kk)=abs(v_eta(kk)-v_eta_an)/abs(v_eta_an);
end
figure;
plot(t,err_csi*r2*1000) % de-normalization absolute error in [m]
figure;
plot(t,err_v_csi)
figure;
plot(t,err_eta*r2*1000) % de-normalization absolute error in [m]
figure;
plot(t,err_v_eta)

figure;
plot(csi,eta)
title("csi-eta HWC")

%
[t,csi,v_csi,eta,v_eta,~,~]=solveHCW([csi0;vplus_csi0;eta0;vplus_eta0;zeta0;vplus_zeta0],linspace(0,Delta_T, 1000));
figure;
axis equal
hold on
mu_Earth = 3.98600418e5;
for kk=1:length(csi)
% f needs to be implemented instead of t(kk)? t_ad*n_ad n_ad=1

[rHt, rH_dott] = get_kepel2hill(r2, 0, t(kk), mu_Earth); % cartesian Hill ref frame, implement f
[rPt, rP_dott] = get_hill2perifocal(rHt,rH_dott,t(kk)); % to perifocal, implement f
[rEt, rE_dott] = get_perifocal2equatorial(rPt,rP_dott,0,0,0); % to equatorial
%plot(rEt(1),rEt(2),"o",Color="b")
% to perifocal implement f
[rP0, rP_dot0] = get_hill2perifocal([csi(kk)*r2+rHt(1);eta(kk)*r2+rHt(2);0],[v_csi(kk)+rH_dott(1);v_eta(kk)+rH_dott(2);0],t(kk));
[rE0, rE_dot0] = get_perifocal2equatorial(rP0,rP_dot0,0,0,0); % to equatorial
plot(rE0(1),rE0(2),"*",Color="b")
% pause(0.1)
end
plot(r2*cos(linspace(0,2*pi,1000)),r2*sin(linspace(0,2*pi,1000)),"r")
hold off%
end
