% this function is the second version of a possible solution

function []=part5_vers2(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)
%% fsolve solution
tspan=[0 Delta_T]; % ToF fixed
coeff_ad=1/2*rho*CD*A/m *r2;
bet=deltaV1;
res=fsolve(@(deltaV_opt) pos_err(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV_opt],[0 Delta_T]),bet);
function err = pos_err(coeff_ad,cond_i,tspan)
    [~,csi,eta,~,~,~,~]= solveHCW_drag(coeff_ad,cond_i,tspan);
    err = [csi(end);eta(end)]' - [0;0];
end
deltaV1=res;
[~,csi,eta,~,v_csi,v_eta,~]=solveHCW_drag(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV1],tspan);
deltaV2=-[v_csi(end);v_eta(end);0];
deltaVtot=norm(deltaV1,2)+norm(deltaV2,2);

fprintf('Final position offset target - chaser with drag relative CCS: %.4f [m]\n', norm([csi(end);eta(end)],2)*r2*1000);
fprintf('using total impulse: %.4f [m/s]\n', deltaVtot*r2*n2*1000);
end