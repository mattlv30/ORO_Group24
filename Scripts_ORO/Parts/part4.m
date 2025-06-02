%% Part 4
% This function solves all the request of the fourth part of the project

function [tau] = part4(L,theta, F)
%% Inputs
% L - link length (vector)
% theta - joints displacement (vector)
% F - contact force on the end effector 

%% Output

%%

% Number od joints
n = length(theta);

% Wrench on the fifth link
w_5 = [0 -F 0 0 0 0]';

% Distance from the centre of the CCS to the base of the fist joint
L0 = 1; %[m]

% Definition fo the vectors connecting each joint to the Centre of mass
% (CoM) of the successive link

g1 = [0 L(1)/2]';
g2 = L(2)/2*[sin(theta(2)) cos(theta(2))]';
g3 = L(3)/2*[sin(theta(2)) cos(theta(2))]';
g4 = L(4)/2*[sin(theta(2)+theta(4)) cos(theta(2)+theta(4))]';
g5 = L(5)/2*[0 1]';

g = [g1 g2 g3 g4 g5];

% Definition of the position of the CoMs of each link in the intertial
% reference frame 

r0 = [0 0]';
r1 = r0 + [0 L0/2]' + g1;
r2 = r1 + L(1)/2*[0 1]' + g2;
r3 = r2 + L(2)/2*[sin(theta(2)) cos(theta(2))]' + g3;
r4 = r3 + L(3)/2*[sin(theta(2)) cos(theta(2))]' + g4;
r5 = r4 + +L(4)/2*[sin(theta(2)+theta(4)) cos(theta(2)+theta(4))]' + g5;

r = [r0 r1 r2 r3 r4 r5];

% Definition of the versors of each revolute joint (note that the versor
% of the joints alligned with the z axis are always constant since we
% assume to be in a 2D space)

e1 = [0 1 0]';
e2 = [0 0 1]';
e3 = [sin(theta(2)) cos(theta(2)) 0]';
e4 = [0 0 1]';
e5 = [sin(theta(2)+theta(4)) cos(theta(2)+theta(4)) 0]';

e = [e1 e2 e3 e4 e5];

figure
axis equal; hold on; grid on;
plot(r(1,:), r(2,:), 'bo')
title('Centre of Mass of each link')
xlabel('x [m]'); ylabel('x [m]');


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

% Rectangualr block daigonal matrix Nd, size 6(n+1)x(6+n)

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

tau = N'*w;

end