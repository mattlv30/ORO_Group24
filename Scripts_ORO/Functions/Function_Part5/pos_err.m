% function used for fsolve optimization in part5

% it is needed to find the zero of this function in order to have a
% "perfect" rendezvous with the perturbated chaser
function err = pos_err(coeff_ad,cond_i,tspan)
        [~,csi,eta,~,~,~,~]= solveHCW_drag(coeff_ad,cond_i,tspan);

        % [0;0] is left just for clarity and it is the final position for a perfect rendezvous
        err = [csi(end);eta(end)]' - [0;0]; 
end