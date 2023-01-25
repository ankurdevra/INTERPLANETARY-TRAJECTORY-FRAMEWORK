function [future_right_ascension,future_decilantion] = future_RA_Dec_relative_rotating_earth_from_orbital_elements(Angular_momentum,eccentricity,inclination,RAAN,Argument_of_perigee,True_anomaly,elapsed_time_hrs)
% The following code calculates the right ascension and declination of an
% earth orbiting satellite FROM ITS INITIAL orbital elements
% after the specified elapsed time.
% It accounts for the oblateness and rotation of earth to calculate the
% final RA nad DEC.
% REQUIRED INPUTS:
% Angular_momentum = angular momentum of orbiting body around earth km^2/sec
% eccentricity = orbital eccentricity
% inclination  = orbital inclination deg
% RAAN = Right ascension of ascending node deg
% Argument_of_perigee = argument of perigee deg
% True_anomaly = true anomaly deg
% elapsed_time_hrs = elapsed time at which we want to calculate the RA and
% DEC hrs
% OUTPUT:
% Right_Ascension = deg, final right ascension after accounting for oblateness and 
% rotation of earth
% Declination = deg, final declination after accounting for oblateness and 
% rotation of earth
%% Creator:- ANKUR DEVRA 
% Develope Date - 2 July 2022
% Iteration 1 -
%% INPUTS
% initial ORBITAL ELEMENTS OF BODY ORBITING EARTH
h = Angular_momentum; % angular momentum of orbiting body around earth km^2/sec
e = eccentricity; %orbital eccentricity
i = inclination; % orbital inclination deg
Big_omega = RAAN; % Right ascension of ascending node deg
Small_omega = Argument_of_perigee; % argument of perigee deg
theta = True_anomaly; % true anomaly deg
elapsed_time_sec = elapsed_time_hrs*3600; % sec, time elapsed

[~,Geocentric_equatorial_state_vectors] = Geocentric_equi_and_perifocal_state_vector_from_orbital_element(h,e,i,Big_omega,Small_omega,theta); % km and km/sec
% geocentric equatorial coords invokes
% Geocentric_equi_and_perifocal_state_vector_from_orbital_element function
% to calculate the initial geocentric equatorial coord from given initial
% orbital elements.

[future_geocentric_equatorial_coordinates,~] = future_geocentric_equatorial_from_initial_accounts_oblateness(Geocentric_equatorial_state_vectors,elapsed_time_hrs); % km and km/sec
% invokes future_geocentric_equatorial_from_initial_accounts_oblateness
% function to calculate the final geocentric equatorial coordinates after
% the defined elapsed time accounting for earth oblateness

omega_earth = (360*(1+(1/365.26)))/(24*3600); % deg/sec angular velocity of earth

theta_earth = omega_earth*elapsed_time_sec; %deg earth rotation in elapsed time

R3_theta_earth = [cosd(theta_earth) sind(theta_earth) 0;-sind(theta_earth) cosd(theta_earth) 0;0 0 1]; % transformation
% DCM from geocentric equatorial XYZ frame to x'y'z' frame used to
% calculate RA and DEC

r_x_bar = R3_theta_earth*[future_geocentric_equatorial_coordinates(1);future_geocentric_equatorial_coordinates(2);future_geocentric_equatorial_coordinates(3)]; % km final geocentric
% equatorial position coordinates after accounting for rotation of earth
% and oblateness.

[Right_Ascension,Declination] = RA_DEC_from_geocentric_equatiorial_position_vector(r_x_bar');% deg
% invokes RA_DEC_from_geocentric_equatiorial_position_vector function to
% calculate the RA and Dec of the final geocentric equatorial position
% coords.

%% OUTPUT
future_right_ascension = Right_Ascension; % deg, final right ascension after accounting for oblateness and 
% rotation of earth
future_decilantion = Declination;% deg, final declination after accounting for oblateness and 
% rotation of earth
end