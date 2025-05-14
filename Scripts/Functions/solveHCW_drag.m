% assuming that the ipothesys are still valid
% V target adimensional Hill CCS Earth [0;1;0]
function [t,csi,eta,zeta,v_csi,v_eta,v_zeta] = solveHCW_drag(coeff_ad,cond_i,tspan)
% standard substitution
% z=[csi,eta,zeta,csi',eta',zeta']

f = @(t,z) [z(4);
    z(5);
    z(6);
    3*z(1)+2*z(5)-coeff_ad*norm([z(4);z(5);z(6)]+[0;1;0]+[-1*z(5);1*z(4);0],2)*(z(4)-1*z(5));
    -2*z(4)-coeff_ad*norm([z(4);z(5);z(6)]+[0;1;0]+[-1*z(5);1*z(4);0],2)*(1+1*z(4)+z(5));
    -z(3)-coeff_ad*norm([z(4);z(5);z(6)]+[0;1;0]+[-1*z(5);1*z(4);0],2)*z(6)];

opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','off');
[t,z]=ode113(f,tspan,cond_i,opts); % check options
csi=z(:,1);
eta=z(:,2);
zeta=z(:,3);
v_csi=z(:,4);
v_eta=z(:,5);
v_zeta=z(:,6);
end