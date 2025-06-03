function [Htransfer, orbit1, orbit2] = plot_Analytical_orbit(a_h,e_h, r1, r2, n1, n2,t_phasing,phi)

span = 500;

% Circular orbits
f = linspace(0, 2*pi, span); 
orbit1 = r1 * [cos(f); sin(f)];
orbit2 = r2 * [cos(f); sin(f)];

% Transfer ellipse
f_h = linspace(pi, 2*pi, span/2);
r_trans = a_h*(1-e_h^2) ./ (1+e_h*cos(f_h));
Htransfer = r_trans .* [cos(f_h) ; sin(f_h)];

% Assumption: Hohmann transfer starts at pi cartesian
f_phasing1=linspace(pi-n1*t_phasing,pi,span); %
f_phasing2=linspace(pi-n2*t_phasing,pi-phi,span); %

figure;
hold on; grid on; axis equal;

% Plot complete circular orbits
plot(orbit1(1,:), orbit1(2,:), 'b', 'LineWidth', 0.5);
plot(orbit2(1,:), orbit2(2,:), 'r', 'LineWidth', 0.5);
plot(0, 0, 'bo', 'MarkerSize', 20, 'MarkerFaceColor','k','HandleVisibility','off');

% Plot phasing orbit 
plot(r1*cos(f_phasing1),r1*sin(f_phasing1),"b", 'LineWidth', 1.5)
plot(r2*cos(f_phasing2),r2*sin(f_phasing2),"r", 'LineWidth', 1.5)
plot(r1*cos(f_phasing1(1)),r1*sin(f_phasing1(1)),"b*", 'LineWidth', 1.5,'HandleVisibility','off')
plot([r2*cos(f_phasing2(1)),r2*cos(f_phasing2(end))],[r2*sin(f_phasing2(1)),r2*sin(f_phasing2(end))],"r*", 'LineWidth', 1.5,'HandleVisibility','off')

% Plot Hohmann transfer
plot(Htransfer(1,:), Htransfer(2,:), 'k', 'LineWidth', 2);
plot(-r1, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor','b','HandleVisibility','off');
plot(r2*cos(2*pi), r2*sin(2*pi), 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r','HandleVisibility','off');
legend('Initial Orbit', 'Target Orbit',"phasing target","phasing chaser", 'Hohmann Transfer');
title('Hohmann transfer manoeuvre');
xlabel('X [km]'); ylabel('Y [km]');

end