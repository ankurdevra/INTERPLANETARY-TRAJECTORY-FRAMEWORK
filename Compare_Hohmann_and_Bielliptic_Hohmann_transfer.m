function [Hohmann_delta_v,Bielliptic_hohmann_delta_v,TOF_hohmann,TOF_bielliptic_hohmann] = Compare_Hohmann_and_Bielliptic_Hohmann_transfer(radius_of_initial_circular_orbit,radius_of_final_circular_orbit,bielliptic_elliptical_trajectory_apogee)
% The following function comapares which orbit transfer type, hohmann or
% bielliptical hohmann, is more efficent to raise initial circular orbit to
% final circular orbit around earth.
% ONLY FOR INITIAL AND FINAL CIRCULAR ORBITS ONLY
% REQUIRED INPUTS:
% radius_of_initial_circular_orbit =  km, initial radius of circular orbit
% bielliptic_elliptical_trajectory_apogee = km, apogee of bielliptic transfer trajectory orbit
% radius_of_final_circular_orbit = km, final radius of desired circular orbit
% OUTPUTS:
% Hohmann_delta_v = km/sec, total delta-v required for hohmann transfer
% Bielliptic_hohmann_delta_v = km/sec, total delta-v required for bielliptic hohmann transfer
% TOF_hohmann = sec, time of flight in hohmann trajectory
% TOF_bielliptic_hohmann = sec, time of flight of body in biellptical transfer trajectory
%% Creator:- ANKUR DEVRA 
% Develope Date - 4 July 2022
% Iteration 1 -
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% INPUTS
r_A = radius_of_initial_circular_orbit; % km, initial radius of circular orbit
r_B = bielliptic_elliptical_trajectory_apogee; % km, apogee of bielliptic transfer trajectory orbit
r_C = radius_of_final_circular_orbit; % km, final radius of desired circular orbit
%% CALCULATIONS
% Bielliptic Hohmann Transfer
% orbit 1 (initial circular orbit)
v_A1 = sqrt(mu_earth/r_A); % km/sec, velocity of body at point A in initial circular orbit 1
% orbit 2 (initial Bielliptic transfer semi-elliptical orbit)
h2 = sqrt(2*mu_earth)*(sqrt((r_A*r_B)/(r_A+r_B))); % km/sec^2; angular momentum of initial Bielliptic transfer semi-elliptical orbit 2
v_A2 = h2/r_A; % km/sec, velocity of body at point A in bielliptical transfer orbit 2
v_B2 = h2/r_B; % km/sec, velocity of body at point B (apogee of bielliptical transfer trajectory) orbit 2
a2 = (1/2)*(r_A + r_B); % km, semimajor axis of initial bielliptic transfer orbit
% orbit 3 (final Bielliptic transfer semi-elliptical orbit)
h3 = sqrt(2*mu_earth)*(sqrt((r_C*r_B)/(r_C+r_B))); % km/sec^2; angular momentum of final Bielliptic transfer semi-elliptical orbit 3
v_B3 = h3/r_B; % km/sec, velocity of body at point B (apogee of bielliptical transfer trajectory) orbit 3
v_C3 = h3/r_C; % km/sec, velocity of body at point C in bielliptic transfer orbit 3
a3 = (1/2)*(r_C+r_B); % km, semimajor axis of final bielliptic transfer orbit
% orbit 4 (final circular orbit)
v_C4 = sqrt(mu_earth/r_C); % km/sec, velocity of body at point C in final circular orbit

total_bielliptic_delta_v = abs(v_A2 - v_A1) + abs(v_B3 - v_B2) + abs(v_C4 - v_C3); % km/sec, total delta-v required for bielliptic hohmann transfer
bielliptic_TOF = (1/2)*(((2*pi)/sqrt(mu_earth))*a2^(3/2) + ((2*pi)/sqrt(mu_earth))*a3^(3/2)); % sec, time of flight of body in biellptical transfer trajectory

% invokes Hohmann_tranfer_from_lower_to_higher_orbit function to calculate
% data for hohmann tranjectory. need to provide spacecraft mass and
% propellant specific impulse of fuel for the following functiont to work.
% arbitrary values of mass and specific impulse will NOT effect resultant
% analysis of this particular code.
[~,~,total_delta_v,~,TOF_hohmann] = Hohmann_tranfer_from_lower_to_higher_orbit(2000,r_A,r_A,r_C,r_C,300);
total_hohmann_delta_v = total_delta_v; % km/sec, total delta-v required for hohmann transfer

if total_bielliptic_delta_v<total_hohmann_delta_v
    X=("Bielliptic Hohmann Transfer more efficient");disp(X);
else
    X=("Hohmann Transfer more efficient");disp(X);
end

%% OUTPUT
Hohmann_delta_v = total_hohmann_delta_v; % km/sec, total delta-v required for hohmann transfer
Bielliptic_hohmann_delta_v = total_bielliptic_delta_v; % km/sec, total delta-v required for bielliptic hohmann transfer
TOF_hohmann = TOF_hohmann; % sec, time of flight in hohmann trajectory
TOF_bielliptic_hohmann = bielliptic_TOF; % sec, time of flight of body in biellptical transfer trajectory
end