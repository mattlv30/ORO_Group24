%% HohmannTransferNumerical
% This function numerically computes the Hohmann transfer. This funciton
% differs from the version 1 becourse it computes also the orbit starting
% from t0 = 0 untile tf = ToF considering both the orbit before and after
% the Hohman transfert. This function uses the Event Function of ODE to
% change integration

% It perform 3 changes of coordinates going from the orbital elements to
% the perifocal reference frame. Then the impulsive
% velocity needed for the Hohmann transfer is added to the orbital
% velocity. Finally, the numerical solution of the problem is computed
% using a numerical solver

%% Inputs
% a, e, i, om, w, f - classical orbital elements (semi major axis,
% eccentricity, inclination, RAAN, argument of periapsis, ture anomaly) at
% the starting point of the transfer
% mu - gravitational constant of the primary 
% ToF - time of flight of the transfer
% h - angular momentum computed at the starting point of the transfer
% DeltaV1 - deltaV to start the manouvre computed with the Hohmann trnasfer
% analitical formula
% DeltaV2 - deltaV to end the manouvre computed with the Hohmann trnasfer
% analitical formula
% t_wait - waiting time for the phasing 

%% Outputs

% t - time evaluation point from the numerical solver (ode113)
% r - vector position in the ECI reference frame [x; y; z]
% r_dot - velocity vector in the ECI reference frame [Vx; Vy; Vz]
% f - ture anomaly

%% Info 

% Authors: Cristian Iacovella, Mattia Li Vigni, Sara Moreira 
% Creation date: 18/05/2025
% Update date:

%% Function 
function [t,r,r_dot,f] = HohmannTransferNumerical_v2(a,e,i, om, w, f, mu, ToF,h, DeltaV1,DeltaV2,  t_wait)

% To cartesian Hill reference frame
[rH0, rH_dot0] = get_kepel2hill(a, e, f, mu); 

% The impulseive manouvre for the hohmann transfer is added here as a deltaV
% The impules is assumed to be tangent to the orbit (theta direction in the
% Hill reference frame)

[rP0, rP_dot0] = get_hill2perifocal(rH0,rH_dot0,f); % to perifocal
[rE0, rE_dot0] = get_perifocal2equatorial(rP0,rP_dot0,om,i,w); % to equatorial

% Definition of the time span (time of flight of the Homann transfer)
orbit_period = 2*pi*sqrt(a^3/mu);

% Definition of the 7 initial condition: 3 position, 3 velocities, 1
% angular momentum 
cond_i = [rE0; rE_dot0; f]';

% Integration of the problem using a function that solves the Kepler
% restricted 2 body problem with a numerical integration (ode113)
KR2BP=@(t,z) [
    z(4);
    z(5);
    z(6);
   -mu*z(1)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
   -mu*z(2)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
   -mu*z(3)/sqrt(z(1)^2+z(2)^2+z(3)^2)^3;
    h/(z(1)^2+z(2)^2+z(3)^2) % in the KR2BP h=cost
    ]; 

%% Waiting orbit

% Trigger function for the first impulse, wait the phasing time to applay
% the first Delta V
wait = @(t,z)  t-t_wait;

wait_event = odeEvent(EventFcn=wait, ...
             Direction="ascending", ...
             Response="stop");

ODE_wait = ode(ODEFcn=KR2BP, ...
    InitialValue=cond_i, ...
    EventDefinition=wait_event, ...
    Solver="ode45", ...
    RelativeTolerance=1e-9, ...
    AbsoluteTolerance=1e-12);

Sol = solve(ODE_wait,0,orbit_period,Refine=100);

% Update the solutions 
t=Sol.Time;
r = Sol.Solution(1:3,:);
r_dot=Sol.Solution(4:6,:);
f = Sol.Solution(end, :);


%% Firts impulse
% Change in the initial condition, sum the first DeltaV

% To add properly the impulse along the tangent direction, consider that
% two parallel vector are proportional 
a = DeltaV1/norm(r_dot(:,end));
DeltaV1_vec = a.*r_dot(:,end);

r_dot_1 = r_dot(:,end)+DeltaV1_vec;

% Update the initial condition after the first impulse 
cond_i1  = [r(:,end); r_dot_1; f(end)]';

% The manouvre after the first impulse last the time of flight 
first_impulse_trigger = @(t,z)  t-(t_wait+ToF);

impulse1_event = odeEvent(EventFcn=first_impulse_trigger, ...
             Direction="ascending", ...
             Response="stop");

% Consider as a starting time the final time of the previous integration
ODE_1 = ode(ODEFcn=KR2BP, ...
    InitialValue=cond_i1, ...
    InitialTime=t(end),...
    EventDefinition=impulse1_event, ...
    Solver="ode45", ...
    RelativeTolerance=1e-9, ...
    AbsoluteTolerance=1e-12);

Sol = solve(ODE_1,t(end),orbit_period,Refine=100);

% Update the solution
t=[t, Sol.Time];
r = [r , Sol.Solution(1:3,:)];
r_dot=[ r_dot, Sol.Solution(4:6,:)];
f = [f, Sol.Solution(end, :)];

%% Second impulse 

% Change in the initial condition, sum the second DeltaV

% Same computation as before, where a minus is added to take into account
% the change due to the CCS
b = abs(DeltaV2)/norm(r_dot(:,end));
DeltaV2_vec = -b.*r_dot(:,end);

r_dot_2 = r_dot(:,end)+DeltaV2_vec;

% Initial condition after the second impulse 
cond_i2  = [r(:,end); r_dot_2; f(end)]';

% This manouvre lasts until the end of one orbital period (convention)
second_impulse_trigger = @(t,z)  t-orbit_period;

impulse2_event = odeEvent(EventFcn=second_impulse_trigger, ...
             Direction="ascending", ...
             Response="stop");

% The initial time is the last instant of the prevoius integration
ODE_2 = ode(ODEFcn=KR2BP, ...
    InitialValue=cond_i2, ...
    InitialTime=t(end),...
    EventDefinition=impulse2_event, ...
    Solver="ode45", ...
    RelativeTolerance=1e-9, ...
    AbsoluteTolerance=1e-12);

Sol = solve(ODE_2,t(end),orbit_period,Refine=100);

% Update the results 
t=[t, Sol.Time];
r = [r , Sol.Solution(1:3,:)];
r_dot=[ r_dot, Sol.Solution(4:6,:)];
f = [f, Sol.Solution(end, :)];

end