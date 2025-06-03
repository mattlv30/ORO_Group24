function [e_x] = e_cross_matrix(Euler_axis)
% This function computes the vectorial product matrix starting from the
% vector Euler_axis

e_x = [0 -Euler_axis(3) Euler_axis(2);
    Euler_axis(3) 0 -Euler_axis(1);
    -Euler_axis(2) Euler_axis(1) 0];

end