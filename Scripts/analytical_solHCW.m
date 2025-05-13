% Description
% input
% output Matrixes
% informations 
% author CI
% creation date 23/03/2025
% update date

function [phi_p, phi_c] = analytical_solHCW(f0)
% z = phi * z0
% phi_p for csi and eta
phi_p=@(f) [4-3*cos(f-f0) sin(f-f0) 0+0*f 2-2*cos(f-f0);
    3*sin(f-f0) cos(f-f0) 0+0*f 2*sin(f-f0);
    -6*(f-f0)+6*sin(f-f0) -2+2*cos(f-f0) 1+0*f -3*(f-f0)+4*sin(f-f0);
    -6+6*cos(f-f0) -2*sin(f-f0) 0+0*f -3+4*cos(f-f0)];
% phi_c for zeta (they aren't coupled)
phi_c=@(f) [cos(f-f0) sin(f-f0);
    -sin(f-f0) cos(f-f0)];
end