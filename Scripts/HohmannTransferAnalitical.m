function [DeltaV1,DeltaV2, DeltaVtot, a_h, e_h, ToF] = HohmannTransferAnalitical(r1, r2, mu)

%HohmannTransfer: This function computes analitically the DeltaV required for a
%Homann trasnfert between two circular orbis of radii r1 and r2, around a
%primary with planetary gravitational constant mu

%The output of the funtion are the two Delta V computed at the perigee and
%apogee of the ellpitical orbit and the time of flight


% Delta V [km/s]
DeltaV1 = sqrt(mu/r1)*(sqrt(2*r2/(r1+r2))-1);
DeltaV2 = sqrt(mu/r2)*(1-sqrt(2*r2/(r1+r2)));
DeltaVtot = abs(DeltaV1) + abs(DeltaV2);

%Semi-major axis of transfer ellipse 
a_h = (r1 + r2)/2; %[km]

%Eccentricity of the transfer ellipse
e_h = (r1-r2)/(r2+r1);
% Time Of Flight
ToF = pi*sqrt((r1+r2)^3/8/mu); %[s]
end