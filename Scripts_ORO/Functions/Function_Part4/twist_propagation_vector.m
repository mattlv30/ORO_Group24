function [p] = twist_propagation_vector(e,g)
% This function computes the twist propagation vector as a function of the
% joint versor e and distance from the following joint centre of mass g,
% according to the formula for the revolute joints

p = [e;
    e_cross_matrix(e)*g];

end