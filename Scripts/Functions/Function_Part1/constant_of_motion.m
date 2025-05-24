function [error_h_rel,error_e_rel,error_E_rel] = constant_of_motion(r,r_dot, t,r1,h_start_H,mu_Earth, e_h, V1 , DeltaV1)
%% Constant of motion :

% 1. Angular momentum h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_vec = cross(r, r_dot);
h_value = vecnorm(h_vec);

% Max relative error between numerical and analitical solution
error_h_rel = max(abs(h_value-ones(1,length(t))*h_start_H))/h_start_H*100;

% 2. Eccentricity vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
e_vec = cross(r_dot, h_vec)/mu_Earth - r./vecnorm(r);
e_value = vecnorm(e_vec);

% Max relative error between numerical and analitical solution
error_e_rel = max(abs(e_value-ones(1,length(t))*e_h))/e_h*100;

% 3. Total energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total orbital energy of the Hohmann transfer (analytical)
E_H = 0.5 * (V1 + DeltaV1)^2 - mu_Earth /r1;

% Total orbital energy of the Hohmann transfer (numerical)
E = 0.5 * vecnorm(r_dot).^2 - mu_Earth ./ vecnorm(r);

% Max relative error between numerical and analitical solution
error_E_rel = max(abs(E - ones(1,length(t))*E_H))/E_H*100;

%% Plot of the results 

figure;
subplot(3,1,1);
plot(t, h_value, 'r', 'LineWidth', 0.5); hold on;
plot(t,ones(1,length(t))*h_start_H, 'b', 'LineWidth', 0.5); grid on;
title('Angular momentum'); ylabel('h $ [\mathrm{km^2}/\mathrm{s}]$', 'Interpreter', 'latex');
lgd_h = legend('Numerical angolar momentum','Analitical angolar momentum');
lgd_h.Position = [0.7 0.79 0.05 0.05]; 
xlabel('Time [s]');

subplot(3,1,2);
plot(t, e_value, 'r', 'LineWidth', 0.5); hold on;
plot(t,ones(1,length(t))*e_h, 'b', 'LineWidth', 0.5); grid on;
title('Eccentricity vector'); ylabel('$e [-]$', 'Interpreter', 'latex');
lgd_e = legend('Numerical eccentricity','Analitical eccentricity');
lgd_e.Position = [0.725 0.5 0.05 0.05]; 
xlabel('Time [s]');

subplot(3,1,3);
plot(t, E, 'r', 'LineWidth', 0.5); hold on;
plot(t,ones(1,length(t))*E_H, 'b', 'LineWidth', 0.5); grid on;
title('Total orbital energy'); ylabel('$\epsilon$', 'Interpreter', 'latex');
lgd_E =legend('Numerical orbital energy','Analitical orbital energy');
lgd_E.Position = [0.725 0.18 0.05 0.05]; 
xlabel('Time [s]');
end