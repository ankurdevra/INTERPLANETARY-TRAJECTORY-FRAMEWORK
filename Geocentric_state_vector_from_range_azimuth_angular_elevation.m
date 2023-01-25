function [Geocentric_state_vector] = Geocentric_state_vector_from_range_azimuth_angular_elevation(range,azimuth,angular_elevation,range_rate,azimuth_rate,angulae_elevation_rate,Local_siderial_time_of_tracking_station_at_time_of_observation,height_above_sea_level_of_tracking_station,latitude_of_tracking_station)
% The following code calculates geocentric position on an orbiting body
% around earth from tracking station observations.
% ALTHOUGH IN THIS CODE ANGLES ARE WITH RESPECT TO DEGREES; RADIAN SHOULD
% BE USED.
% CODE NEEDS REFINEMENT AS CONFUSION IN DEG AND RAD
% INCOSISTENT VALUES OF GEOCENTRIC VELOCITY WITH CHOICE OF ANGLE UNIT, DEG
% OR RAD
% FUTURE RATIFICATION REQUIRED
% REQUIRED INPUTS:
% range; % km slant range (distance) distance of satellite or celestial body (B) from observer O (tracking site)
% azimuth; % azimuth, deg, measured clockwise from due north (0<=A<=360) deg
% angular_elevation; % deg, angular elevation or altitude measure from the horizontal to the line of sight on the body B
% (-90<=a<=90) deg
% range_rate; % deg/sec, rate of change of range wrt time
% azimuth_rate; % deg/sec, rate of change of azimuth wrt time
% angulae_elevation_rate; % deg/sec, rate of change of angular elevation wrt time
% Local_siderial_time_of_tracking_station_at_time_of_observation; % deg, local siderial time of tracking station at the time of observation
% height_above_sea_level_of_tracking_station; % km, height of tracking station above sea level
% latitude_of_tracking_station; % deg, latitude of tracking station
% OUTPUT
% Geocentric_state_vector = km and km/sec, geocentric state vector of body
% B from tracking station measurements
%% Creator:- ANKUR DEVRA 
% Develope Date - 4 July 2022
% Iteration 1 -
%% Starting data
R_earth = 6378.137; % equatorial radius of earth km
f = 0.003353; % flattening factor f of earth
omega_earth = (360*(1+(1/365.26)))/(24*3600); % deg/sec angular velocity of earth
%% INPUTS
% Topocentric Horizon system
% Quantities wrt tracking station on earth
rho = range; % km slant range (distance) distance of satellite or celestial body (B) from observer O (tracking site)
A = azimuth; % azimuth, deg, measured clockwise from due north (0<=A<=360) deg
a = angular_elevation; % deg, angular elevation or altitude measure from the horizontal to the line of sight on the body B
% (-90<=a<=90) deg
rho_dot = range_rate; % deg/sec, rate of change of range wrt time
A_dot = azimuth_rate; % deg/sec, rate of change of azimuth wrt time
a_dot = angulae_elevation_rate; % deg/sec, rate of change of angular elevation wrt time
theta = Local_siderial_time_of_tracking_station_at_time_of_observation; % deg, local siderial time of tracking station at the time of observation
H = height_above_sea_level_of_tracking_station; % km, height of tracking station above sea level
phi = latitude_of_tracking_station; % deg, latitude of tracking station
%% CALCULATION
R1_and_R2 = ((((R_earth/(sqrt(1-(2*f - f^2)*(sind(phi))^2)))) + H)*cosd(phi)).*[cosd(theta);sind(theta);0]';
R3 = ((((((R_earth)*(1-f^2))/(sqrt(1-(2*f - f^2)*(sind(phi))^2)))) + H).*sind(phi)).*[0;0;1]';
R = (R1_and_R2 + R3)';
%R = ((((R_earth/(sqrt(1-(2*f - f^2)*(sind(phi))^2)))) + H)*cosd(phi)).*[cosd(theta);sind(theta)]' + ((((((R_earth)*(1-f^2))/(sqrt(1-(2*f - f^2)*(sind(phi))^2)))) + H).*sind(phi))% km, Geocentric position vector of observer (tracking station)
% directed from center of earth to tracking station
delta = asind(((cosd(phi))*(cosd(A))*(cosd(a))) + ((sind(phi))*(sind(a)))); % deg, topocentric delination, declination of body (B) wrt tracking station
h = acosd((((cosd(phi))*(sind(a))) - ((sind(phi))*(cosd(A))*(cosd(a))))/(cosd(delta))); %deg, hour angle (h) is the
% angular distance between the object and the local meridian. IF h is
% positive, the object is west of meridian, if h is negative the object is
% east of meridian
% resolving quadrant ambiguity
if 0<A && A<180
    h = 360-h;
elseif 180<=A && A<=360
    h = h;
end
alpha = theta-h; % deg, topocentric right ascension
rho_hat = cosd(delta).*[cosd(alpha);sind(alpha);0]' + sind(delta).*[0;0;1]'; % direction cosine unit vector. Unit vector of body B wrt topocentric horizon coord system in line of sight direction
r = (R + rho.*rho_hat')'; % km, geocentric position vector of body B
Big_omega = [0;0;omega_earth]'; % deg/sec angular velocity of earth vector
R_dot = cross(Big_omega,R); %km/sec, inertial velocity of tracking station (observer) 
delta_dot = ((1/(cosd(delta))))*((-A_dot*(cosd(phi))*(sind(A))*(cosd(a))) + (a_dot*((((sind(phi))*(cosd(a)))-(((cosd(phi))*(cosd(A))*(sind(alpha)))))))); % deg/sec, declination rate
alpha_dot = omega_earth + ((A_dot*(cosd(A))*(cosd(a))) - (a_dot*(sind(A))*(sind(a))) + (delta_dot*(sind(A))*(cosd(a))*(tand(delta))))/(((cosd(phi))*(sind(a))) - ((sind(phi))*(cosd(A))*(cosd(a)))); %deg/sec, right ascension rate
rho_hat_dot = [(-alpha_dot*(sind(alpha))*(cosd(delta))) - (delta_dot*(cosd(alpha))*(sind(delta)));(alpha_dot*(cosd(alpha))*(cosd(delta))) - (delta_dot*(sind(alpha))*(sind(delta)));delta_dot*cosd(delta)]';% deg/sec direction cosine rate vector
v = (R_dot.*[1;1;1]' + (rho.*(rho_hat_dot)).*[1;1;1]' + (rho_dot.*(rho_hat)).*[1;1;1]'); % km/sec, geocentric velocity vector of body B
%% OUTPUT
Geocentric_state_vector = [r,deg2rad(v)]; % km and km/sec, geocentric state vector of body B
end