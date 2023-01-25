function [Angular_momentum,Inclination,Eccentricity,RAAN,Argument_of_Perigee,True_anomaly,Radial_velocity,Time_period_of_orbit_hrs,Radius_perigee,Radius_apogee,Semimajor_axis,Semilatus_rectum,Eccentric_anomaly,Mean_anomaly,Orbits_in_a_day] = Orbital_element_from_triad_of_geocentric_position_vectors_Gibbs(r1_vec,r2_vec,r3_vec)
% The following function calculates orbital elements of an body orbiting
% earth whose three geocentric position vectors are available.
% It used Gibbs method of preliminary orbit determination
% The anomalies are with regards to the second [r2_vec,v2_vec] geocentric
% state vector
% REQUIRED INPUTS:
% r1_vec = [1x3] GEOCENTRIC position vector of first observation
% r2_vec = [1x3] GEOCENTRIC position vector of second observation
% r3_vec = [1x3] GEOCENTRIC position vector of third observation
% OUTPUTS:
% Angular_momentum = specific angular momentum km^2/sec
% Inclination = inclination of orbit is degrees, if i>90 retrograde orbit
% Eccentricity =  eccentricity
% RAAN =  Right Ascension of Ascending Node, degrees
% Argument_of_Perigee =  argument of perigee degrees
% True_anomaly =  true anomaly degrees
% Radial_velocity =  km/sec (v_r) if v_r>0 trajectory away from preigee if v_r<0 towards perigee
% Time_period_of_orbit_hrs = Time period of orbit in hrs
% Radius_perigee =  radius of perigee in km
% Radius_apogee =  radius of apogee in km
% Semimajor_axis =  semi-major axis in km
% Semilatus_rectum =  semi-latus rectum in km
% Eccentric_anomaly =  Eccentric anomaly in degrees
% Mean_anomaly =  Mean anomaly in degrees
% Orbits_in_a_day =  number of orbits in a day 
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% Creator:- ANKUR DEVRA 
% Develope Date - 2 July 2022
% Iteration 1 -
%% CALCULATIONS
r1_mag = norm(r1_vec); % km, maginute of first position vector
r2_mag = norm(r2_vec); % km, maginute of second position vector
r3_mag = norm(r3_vec); % km, maginute of third position vector

C12_vec = cross(r1_vec,r2_vec); %km^2
C23_vec = cross(r2_vec,r3_vec); %km^2
C31_vec = cross(r3_vec,r1_vec); %km^2

N_vec = r1_mag.*C23_vec + r2_mag.*C31_vec + r3_mag.*C12_vec; % km^3
N_mag = norm(N_vec); % km^3

D_vec = C12_vec+C23_vec+C31_vec; % km^2
D_mag = norm(D_vec); % km^2

S_vec = r1_vec.*(r2_mag-r3_mag) + r2_vec.*(r3_mag-r1_mag) + r3_vec.*(r1_mag-r2_mag); % km^2

v2_vec = (sqrt((mu_earth)/(N_mag*D_mag)))*(((cross(D_vec,r2_vec))/r2_mag)+S_vec); % km/sec velocity vector v2 corresponding to r2 position

geocentric_state_vector = [r2_vec,v2_vec]; % km and km/sec geocentric state vector of orbiting body
%% OUTPUT
% invokes Orbital_elements_from_State_vectors function to calulate the
% orbital elements with respect to second observation
[Angular_momentum,Inclination,Eccentricity,RAAN,Argument_of_Perigee,True_anomaly,Radial_velocity,Time_period_of_orbit_hrs,Radius_perigee,Radius_apogee,Semimajor_axis,Semilatus_rectum,Eccentric_anomaly,Mean_anomaly,Orbits_in_a_day] = Orbital_elements_from_State_vectors(geocentric_state_vector);
end