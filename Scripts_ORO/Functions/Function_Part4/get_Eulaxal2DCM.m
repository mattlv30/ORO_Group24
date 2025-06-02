function [DCM] = get_Eulaxal2DCM(e,alpha)
% This function computes the rotation matrix associated to the vector e and
% the angle of rotation alpha

DCM = cos(alpha)*eye(3) + (1-cos(alpha)).*e*e' - sin(alpha)*e_cross_matrix(e);

end