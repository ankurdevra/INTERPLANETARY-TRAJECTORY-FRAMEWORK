function [Delta_v,Time_to_Change_Plane] = Plane_Change_Maneuver(Current_Apogee,Current_Perigee,Current_Inclination,Desired_Apogee,Desired_Perigee,Desired_Inclination,True_Anomaly_at_Initiation)
% The following function determines the delta v required and time it takes
% for plane change maneuver from one orbit to another to happen.
% REQUIRED INPUTS:
% Current_Apogee = Current apogee of spacecraft in km
% Current_Perigee = Current perigee of spacecraft in km
% Current_Inclination = Current inclination of spacecraft in deg
% Desired_Apogee = Desired apogee of spacecraft in km
% Desired_Perigee = Desired perigee of spacecraft in km
% Desired_Inclination = Desired inclination of spacecraft in deg
% True_Anomaly_at_Initiation = Ture anomaly at which we want to initiate
% the plane change maneuver deg
% OUTPUTS:
% Delta_v = change in velocity delta v required for plane change manuever
% km/sec
% Time_to_Change_Plane = Time it takes for plane change to happen sec
%% Creator:- ANKUR DEVRA 
% Develope Date - 24 March 2022
% Iteration 1 -
%% Starting data
R_earth = 6378.137; % km Earth Equatorial Radius
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
% current
ra1 = Current_Apogee; % km current orbit apogee
rp1 = Current_Perigee; % km current orbit perigee
i1 = Current_Inclination; % deg current orbit inclination
% desired
ra2 = Desired_Apogee; % km desired orbit apogee
rp2 = Desired_Perigee; % km desired orbit perigee
i2 = Desired_Inclination; % deg desired orbit inclination
theta = True_Anomaly_at_Initiation; % deg true anomaly on which to exicute plane change
delta = i2-i1; % change in inclination
%% Current orbit calculations 
% Note:- calls a function with plots, do not run it in loop until
% commenting out the plot.
[Eccentricity_Current,Angular_momentum_Current] = Orbital_Data_From_Apogee_and_Perigee((-R_earth+rp1),(-R_earth+ra1));
e1 = Eccentricity_Current; % eccentricity of current orbit
h1 = Angular_momentum_Current; % km^2/sec angular momentum of current orbit
vt1 = ((mu_earth)/(h1))*(1+(e1*cosd(theta))); % km/sec azimuthal/transverse component of current orbit velocity at given true anomaly
vr1 = ((mu_earth)/(h1))*((e1*sind(theta))); % km/sec radial component of velocity current orbit at given true anomaly
r1 = ((h1)^2/(mu_earth))*(1/(1+(e1*cosd(theta)))); % km radial position of manuevering point for current orbit

%% Desired orbit calculations 
[Eccentricity_Desired,Angular_momentum_Desired] = Orbital_Data_From_Apogee_and_Perigee((-R_earth+rp2),(-R_earth+ra2));
e2 = Eccentricity_Desired; % eccentricity of desired orbit
h2 = Angular_momentum_Desired; % km^2/sec angular momentum of desired orbit
vt2 = ((mu_earth)/(h2))*(1+(e2*cosd(theta))); % km/sec azimuthal/transverse component of desired orbit velocity at given true anomaly
vr2 = ((mu_earth)/(h2))*((e2*sind(theta))); % km/sec radial component of desired orbit velocity at given true anomaly
r2 = ((h2)^2/(mu_earth))*(1/(1+(e2*cosd(theta)))); % km radial position of manuevering point for desired orbit
%% Output
Delta_v = sqrt((vr2-vr1)^2 + (vt1)^2 + (vt2)^2 - (2*vt1*vt2*cosd(delta))); % km/sec required delta v for plane change manuver
Time_to_Change_Plane = (abs(r2-r1))/Delta_v; % sec time it takes to change the plane of orbit
%%
% %% Desired orbit calculations 
% [Eccentricity_Desired,Angular_momentum_Desired,velo3] = Orbital_Data_From_Apogee_and_Perigee((-R_earth+rp1),(-R_earth+ra2));
% velo3
% e2_ = Eccentricity_Desired; % eccentricity of desired orbit
% h2_ = Angular_momentum_Desired; % km^2/sec angular momentum of desired orbit
% vt2_ = ((mu_earth)/(h2))*(1+(e2*cosd(theta))); % km/sec azimuthal/transverse component of desired orbit velocity at given true anomaly
% vr2_ = ((mu_earth)/(h2))*((e2*sind(theta))); % km/sec radial component of desired orbit velocity at given true anomaly
% r2_ = ((h2)^2/(mu_earth))*(1/(1+(e2*cosd(theta)))); % km radial position of manuevering point for desired orbit
end