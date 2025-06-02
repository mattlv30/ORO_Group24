%% Function used for fsolve optimization in part5


%% Input 

% coeff_ad - coefficient takeing into account all the constant parameters
% for the computation of the drag 
% cond_i - initial condition
% tspan - timespan

% 
% it is needed to find the zero of this function in order to have a
% "perfect" rendezvous with the perturbated chaser
function err = pos_err(coeff_ad,cond_i,tspan)
        [~,csi,eta,~,~,~,~]= solveHCW_drag(coeff_ad,cond_i,tspan);

        % [0;0] is left just for clarity and it is the final position for a perfect rendezvous
        err = [csi(end);eta(end)]' - [0;0]; 
end