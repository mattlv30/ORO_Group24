%% Hill reference frame to perifocal reference frame 
% Description state vector from Hill frame to perifocal frame

%% Input 
% rH rH_dot Hill ref frame state vector, true anomaly f

%% Output 
% rP rP_dot perifocal frame state vector

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 09/05/2025
% Update date:

function [rP,rP_dot] = get_hill2perifocal(rH,rH_dot,f)
C_PH = [cos(f) -sin(f) 0; sin(f) cos(f) 0; 0 0 1];
rP=C_PH*rH;
rP_dot=C_PH*rH_dot;
end