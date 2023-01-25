function [Future_State_Vector,Eccentricity] = Future_State_Vector_From_Initial_Using_Lagrange_Coefficients(Initial_State_Vector,Change_in_True_Anomaly)
% Calculates future state vector of an earth orbitting body from its
% initial state vectors and change in true anomaly ( in degrees ) using
% lagrange coefficent method. Also finds the eccentricity of resultant
% trajectory.
% REQUIRED INPUTS:
% Initial_State_Vector = row vector containing the initial state vector of
% body in km and km/sec
% Change_in_True_Anomaly = the change in true anomaly from two subsequent
% measurements, in degrees
% change in true anomaly (degrees) = initial true anomaly - future true anomaly(future state vector is found at this true anomaly) 
% OUTPUTS:
% Future_State_Vector = final state vector km and km/sec
% Eccentricity = eccentricity of trajectory
% Note- all angle measurements are in degrees
%% Creator:- ANKUR DEVRA 
% Develope Date - 23 March 2022
% Iteration 1 -
%% Starting data
format long;
mu_earth = 398600.4418;% Gravitational parameter of earth in km^3/s^2
delta_theta = Change_in_True_Anomaly; % given change in true anomaly in degrees
r0 = Initial_State_Vector(1:3);% initial position vector in km
v0 = Initial_State_Vector(4:6);% initial velocity vector in km/sec
%% Calculations
r0_mag = norm(r0); % initial distance km
v0_mag = norm(v0); % initial velocity km/sec
vr0_mag = (dot(v0,r0))/r0_mag; % radial component of initial velocity km/sec
h = r0_mag*sqrt((v0_mag)^2 - (vr0_mag))^2; % angular momentum km^2/sec
r = (h^2/mu_earth) * (1/(1 + (((h^2/(mu_earth*r0_mag))-1)*cosd(delta_theta))-((h*vr0_mag)/mu_earth)*sind(delta_theta))); % final distance of body km
f = 1-((mu_earth*r)/h^2)*(1-cosd(delta_theta)); % lagrange coefficient dimensionless
g = ((r*r0_mag)/h)*sind(delta_theta); % lagrange coefficient dimensionless
f_dot = ((mu_earth*(1-cosd(delta_theta)))/(h*sind(delta_theta)))*((((mu_earth/h^2)*(1-cosd(delta_theta))-(1/r0_mag)-(1/r)))); % lagrange coefficient dimensionless
g_dot = (1+(f_dot*g))/f; % lagrange coefficient dimensionless
r_future = f.*(r0) + g.*(v0); % final position vector km
v_future = f_dot.*(r0) + g_dot.*(v0); % final velocity vector km/sec
e =sqrt(((h^2/(mu_earth*r0_mag))-1)^2 + (vr0_mag*h/mu_earth)^2); % eccentricity of trajectory
%% Outputs
Future_State_Vector = [r_future,v_future]; % final state row (1X6) vector matrix km and km/sec
Eccentricity = e; % eccentricity of trajectory
end
