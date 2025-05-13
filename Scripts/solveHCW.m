% Description integrates the normalized Hill-Clohessy-Wiltshire equations
% input tspan, intial condition [cis0,v_csi0,eta0,v_eta0,zeta0,v_zeta0]
% output 
% informations 
% author CI
% creation date 23/03/2025
% update date

function [t,csi,v_csi,eta,v_eta,zeta,v_zeta] = solveHCW(cond_i,tspan)
% standard substitution
% z=[csi,csi',eta,eta',zeta,zeta']
f = @(t,z) [z(2);
    3*z(1)+2*z(4);
    z(4);
    -2*z(2);
    z(6);
    -z(5)];
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');
[t,z]=ode113(f,tspan,cond_i,opts); % check options
csi=z(:,1);
v_csi=z(:,2);
eta=z(:,3);
v_eta=z(:,4);
zeta=z(:,5);
v_zeta=z(:,6);
end