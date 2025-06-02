% this function is the second version of a possible solution
% ToF is fixed (using the one calculated from the unperturbated rendezvous)

function deltaVtot=part5_vers2(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)



%% fsolve solution
tspan=[0 Delta_T]; % ToF fixed
coeff_ad=1/2*rho*CD*A/m *r2;
bet=deltaV1; % ideal deltaV1 is chosen as first bet

% hide fsolve output in command window
options = optimoptions("fmincon", "Display", "off");
warning("off", "all");

% find 1st impulse such that zero of the position at the end of the rendezvous
res=fsolve(@(deltaV_opt) pos_err(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV_opt],[0 Delta_T]),bet,options);
warning("on","all")

deltaV1=res; % to try the calculated first impulse
[~,csi,eta,~,v_csi,v_eta,~]=solveHCW_drag(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV1],tspan);
deltaV2=-[v_csi(end);v_eta(end);0]; % second impulse
deltaVtot=norm(deltaV1,2)+norm(deltaV2,2); % total impulse

% disp(" ==================== Part 5 - Optimisation of Delta V and ToF ==================== ")
% disp(" ")
% 
% fprintf('Final position offset target - chaser with drag relative CCS: %.4f [m]\n', norm([csi(end);eta(end)],2)*r2*1000);
% fprintf('ToF: %.4f [s]\n', Delta_T/n2);
% fprintf('using total impulse: %.4f [m/s]\n', deltaVtot*r2*n2*1000);
% 
% disp(" ")
end