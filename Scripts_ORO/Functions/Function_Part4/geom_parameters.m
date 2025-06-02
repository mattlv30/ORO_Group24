function [r,g, e,j, links] = geom_parameters(L,theta, L0)

% This function defines the geometrical parameter of the robotic
% manipulator

%% Input 
% L - lenght of the links
% theta - angular displacement of the joints

%% Output 
% r - position vector of the centres of mass (CoMs) of the links in the
% inertial reference frame 
% g - vector from the centre of the joint to the CoM of the successive link
% e - versor of the joints 
% j - position vector of the joints in the inertial reference frame 




% Definition fo the vectors connecting each joint to the Centre of mass
% (CoM) of the successive link

g1 = [0 L(1)/2]';
g2 = L(2)/2*[sin(theta(2)) cos(theta(2))]';
g3 = L(3)/2*[sin(theta(2)) cos(theta(2))]';
g4 = L(4)/2*[sin(theta(2)+theta(4)) cos(theta(2)+theta(4))]';
g5 = L(5)/2*[sin(theta(2)+theta(4)) cos(theta(2)+theta(4))]';

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

% Position of the joints in the inertial reference frame 

j = r(:, 2:end) - g;

% Links in the inertial reference frame 

links = [j j(:,end)+L(4)*[sin(theta(2)+theta(4)) cos(theta(2)+theta(4))]'];



end