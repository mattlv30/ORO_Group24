function [] = plot_numerical(Htransfer,r, a_h,e_h,f)

%% Plot of the orbits

figure
axis equal; hold on; grid on;
plot(0, 0, 'bo', 'MarkerSize', 20, 'MarkerFaceColor','k','HandleVisibility','off');

% Plot the analitical Hohmann transfer
plot(Htransfer(1,:), Htransfer(2,:), 'r-', 'LineWidth', 2);

% Plot the numerical Hohmann transfer
plot(r(1,:),r(2,:), 'b-', 'LineWidth', 1);
title('Numerical and analitical Hohmann transfer')
legend("Numerical transfer ","Analitical transfer");
xlabel('X [km]'); ylabel('Y [km]');

%% Plot of relative error

% Compute the relative error between the analitical and numerical orbit
r_trans_analitic = a_h*(1-e_h^2) ./ (1+e_h*cos(f));

figure; grid on;
error_rel = (vecnorm(r)-r_trans_analitic)./r_trans_analitic;
plot(f, error_rel)
title('Relative error')
xlabel('True anomaly [rad]'); ylabel('Relative error [%]');

end