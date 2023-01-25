function [Perifocal_state_vectors,Geocentric_equatorial_state_vectors] = Geocentric_equi_and_perifocal_state_vector_from_orbital_element(Angular_momentum,eccentricity,inclination,RAAN,Argument_of_perigee,True_anomaly)
% The following code converts orbital elements of a body around earth to perifocal and geocentric
% equatorial frame coordinates.
% REQUIRED INPUTS:
% Angular_momentum = angular momentum of orbiting body around earth km^2/sec
% eccentricity = orbital eccentricity
% inclination  = orbital inclination deg
% RAAN = Right ascension of ascending node deg
% Argument_of_perigee = argument of perigee deg
% True_anomaly = true anomaly deg
% OUTPUT:
% Perifocal_state_vectors = [1X6] State vector in perifocal frame km and km/sec
% Geocentric_equatorial_state_vectors = [1X6] State vector in geocentric equatorial frame km and km/sec
%% Creator:- ANKUR DEVRA 
% Develope Date - 2 July 2022
% Iteration 1 -
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% INPUTS
% ORBITAL ELEMENTS OF BODY ORBITING EARTH
h = Angular_momentum; % angular momentum of orbiting body around earth km^2/sec
e = eccentricity; %orbital eccentricity
i = inclination; % orbital inclination deg
Big_omega = RAAN; % Right ascension of ascending node deg
Small_omega = Argument_of_perigee; % argument of perigee deg
theta = True_anomaly; % true anomaly deg
%% CALCULATIONS
Perifocal_position_vector = (((h^2)/(mu_earth))*(1/(1+e*cosd(theta)))).*[cosd(theta);sind(theta);0]; % km [3X1] position vector in perifocal frame
Perifocal_velocity_vector = ((((mu_earth)/h))).*[-sind(theta);e+cosd(theta);0]; % km/sec [3X1] velocity vector in perifocal frame

Q_X_x_bar = [-sind(Big_omega)*cosd(i)*sind(Small_omega)+cosd(Big_omega)*cosd(Small_omega) cosd(Big_omega)*cosd(i)*sind(Small_omega)+sind(Big_omega)*cosd(Small_omega) sind(i)*sind(Small_omega);
             -sind(Big_omega)*cosd(i)*cosd(Small_omega)-cosd(Big_omega)*sind(Small_omega) cosd(Big_omega)*cosd(i)*cosd(Small_omega)-sind(Big_omega)*sind(Small_omega) sind(i)*cosd(Small_omega);
             sind(Big_omega)*sind(i) -cos(Big_omega)*sind(i) cosd(i)]; % DCM to transform from geocentric equatiorial frame to perifocal frame

Q_x_bar_X = Q_X_x_bar'; % DCM to transform from perifocal to geocentric equatorial frame

Geocentric_equatorial_position_vector = Q_x_bar_X*Perifocal_position_vector; %km [3X1] position vector in geocentric equatorial frame
Geocentric_equatorial_velocity_vector = Q_x_bar_X*Perifocal_velocity_vector; %km/sec [3X1] velocity vector in geocentric equatorial frame
%% OUTPUT
Perifocal_state_vectors = [Perifocal_position_vector',Perifocal_velocity_vector']; % km and km/sec [1X6] state vector in perifocal frame
Geocentric_equatorial_state_vectors = [Geocentric_equatorial_position_vector',Geocentric_equatorial_velocity_vector']; % km and km/sec [1X6] state vector in geocentric equatorial frame
end