%% Perfiocal to equatorial 
% Description from state vector in perifocal frame to equatorial frame
% input rP rP_dot perifocal, om right ascension of the ascending node OMEGA
%       w argument of periapsis omega, i inclination i
% output rE rP_dot equatorial frame state vector
% informations elementary rotation angle f about angular momentum direction
%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 09/05/2025
% Update date:

function [rE,rE_dot] = get_perifocal2equatorial(rP,rP_dot,om,i,w)

C_EP=[cos(om)*cos(w)-sin(om)*cos(i)*sin(w) -cos(om)*sin(w)-sin(om)*cos(i)*cos(w) sin(om)*sin(i);
     sin(om)*cos(w)+cos(om)*cos(i)*sin(w) -sin(om)*sin(w)+cos(om)*cos(i)*cos(w) -cos(om)*sin(i);
     sin(i)*sin(w) sin(i)*cos(w) cos(i)];

rE=C_EP*rP;
rE_dot=C_EP*rP_dot;
end