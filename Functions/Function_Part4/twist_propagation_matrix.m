function [Bij] = twist_propagation_matrix(c_ij)
% This function creates the twist propagation matrix Bij starting from the
% vector connecting the centers of mass of two successive links

Bij = [eye(3) zeros(3);
    e_cross_matrix(c_ij) eye(3)];

end