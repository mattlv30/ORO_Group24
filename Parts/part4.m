%% Part 4
% This function solves all the request of the fourth part of the project

function [tau] = part4(L,theta, w_5)
%% Inputs
% L - link length (vector)
% theta - joints displacement (vector)
% F - contact force on the end effector 

%% Output

% tau - generalised force matrix representing the combined base and
% manipulator generalised forces (first 6 components are the resultant
% moment and forces on the base in inertial coordinates, while the other
% are the actuating force or torques on the i-th joint of the manipulator)

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 30/05/2025
% Update date: 02/06/2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definition of geometrical parameters

% Number od joints
n = length(theta);

% Distance from the centre of the CCS to the base of the fist joint
L0 = 1; %[m]

[r,g,e,j_pos, links] = geom_parameters(L,theta, L0);


%% System Jacobian Matrix

% Square block lower triangular matrix Nl, size 6(n+1)x6(n+1)

Nl = zeros(6*(n+1),6*(n+1));

for i = 1:n+1
    for j=1:n+1
        if i==j
            Bii = eye(6);
            Nl(6*i-5:6*i, 6*j-5:6*j) = Bii;
        elseif i>j
            c_ij = r(:,j) - r(:,i);
            Bij = twist_propagation_matrix([c_ij;0]);
            Nl(6*i-5:6*i, 6*j-5:6*j) = Bij;
        end
    end
end

% Rectangualr block diagonal matrix Nd, size 6(n+1)x(6+n)

Nd = zeros(6*(n+1),(6+n));
P0 = eye(6);

Nd(1:6, 1:6) = P0;


for i = 1:n
    pi = twist_propagation_vector(e(:,i), [g(:,i) ; 0]);
    Nd(6+6*i-5:6+6*i,6+i) = pi;
end

% System Jacobian Matrix N

N = Nl*Nd;

% Wrench of the manipulator 
w = zeros(6*(n+1),1);
w (end-5:end,:) = w_5;

%% Computation of the generalised force matrix

tau = N'*w;

%% Plot and display 

figure
hold on; grid on; %axis equal;
plot(links(1,:), links(2,:), 'k', 'LineWidth', 2)
plot(r(1,2:end), r(2,2:end), 'bo', 'LineWidth', 1.5)
plot(j_pos(1,:), j_pos(2,:), 'r*', 'LineWidth', 1.5)
plot([-L0/2 -L0/2 L0/2 L0/2], [0 L0/2 L0/2 0], 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
plot(links(1,end), links(2,end), 'go', 'MarkerSize', 10, 'MarkerFaceColor','b','HandleVisibility','off');
legend('Links', 'CoM of Links', 'Joints', 'Base')

title('Robotic manipulator configuration')
xlabel('x [m]'); ylabel('y [m]');
xlim([-1 1]); ylim([0 1.6]);

disp(" ==================== Part 4 - Robotic Manipulator ==================== ")
disp(" ")

fprintf('Torque on the base (along z axis): %.4f Nm \n', tau(3));
fprintf('Force on the base (along x axes): %.4f N \n', tau(4));
fprintf('Force on the base (along x axes): %.4f N \n', tau(5));
disp(" ")
fprintf('Torque on the 2nd joints: %.4f Nm \n',  tau(8));
fprintf('Torque on the 4th joints: %.4f Nm \n',  tau(10));
disp(" ")


end