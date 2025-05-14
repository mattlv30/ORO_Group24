% Description state vector from Hill frame to perifocal frame
% input rH rH_dot Hill ref frame state vector, true anomaly f
% output rP rP_dot perifocal frame state vector
% informations elementary rotation angle f about angular momentum direction
% author CI
% creation date 18/03/2025
% update date

function [rP,rP_dot] = get_hill2perifocal(rH,rH_dot,f)
C_PH = [cos(f) -sin(f) 0; sin(f) cos(f) 0; 0 0 1];
rP=C_PH*rH;
rP_dot=C_PH*rH_dot;
end