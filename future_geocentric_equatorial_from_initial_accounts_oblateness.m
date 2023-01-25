function [future_geocentric_equatorial_coordinates,future_perifocal_coordinates] = future_geocentric_equatorial_from_initial_accounts_oblateness(initial_geocentric_equatorial_coordinates,delta_t_hrs) 
% The following code calculates the future geocentric equatorial and perifocal coordinates after specified elapsed time
% from initial geocentric equatorial coordinates and accounts for the J2
% perturbation (oblatenss) of earth.
% REQUIRED INPUTS:
% initial_geocentric_equatorial_coordinates = [1X6] initial State vector in geocentric equatorial frame km and km/sec
% delta_t_hrs = elapsed time in hrs
% OUTPUT:
% future_geocentric_equatorial_coordinates = km and km/sec [1X6] future geocentric equatorial coordinates account for earth oblateness
% future_perifocal_coordinates = km and km/sec [1X6] future perifocal coordinates account for earth oblateness
%% Creator:- ANKUR DEVRA 
% Develope Date - 2 July 2022
% Iteration 1 -
%% Starting data
format long;
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
J2_earth = 1.08263*10^(-3); % J2, second zonal harmonics of earth, used to account for oblateness of earth
R_earth = 6378.137; % equatorial radius of earth km
%% INPUTS
state_vector = initial_geocentric_equatorial_coordinates; % km and km/sec, initial state vector of earth orbiting body
delta_t = 3600*delta_t_hrs; % elapsed time in sec
%% CALCULATIONS
[Angular_momentum,Inclination,Eccentricity,RAAN,Argument_of_Perigee,~,~,Time_period_of_orbit_hrs,~,~,Semimajor_axis,~,Eccentric_anomaly,~,~] = Orbital_elements_from_State_vectors(state_vector); % invokes Orbital_elements_from_State_vectors function
% to find orbital elements from state vector
h = Angular_momentum; % specific angular momentum km^2/sec
i = Inclination; % inclination of orbit is degrees, if i>90 retrograde orbit
e = Eccentricity; % eccentricity
Big_omega = RAAN; % initial Right Ascension of Ascending Node, degrees
Small_omega = Argument_of_Perigee; % initial argument of perigee degrees
%theta = True_anomaly; % initial true anomaly degrees
a = Semimajor_axis; % semi-major axis in km
Time_period_sec = Time_period_of_orbit_hrs*3600; % Time period of orbit in sec
n = (2*pi)/Time_period_sec; % rad/sec mean motion
E = (Eccentric_anomaly)*(pi/180); % initial eccentric anomaly rad
% Solving keplers equation to calculate time since perigee passage at
% inital epoch
t0 = (E-e*sin(E))/n;% sec time since perigee passage at inital epoch
tf = t0+delta_t; % sec final time after elapsed time
np = tf/Time_period_sec; % number of period since passing perigee in first orbit;that is number of orbits
% completed after elpased time since first perigee passage.
tp = (np - floor(np))*Time_period_sec; % sec time since perigee passage in the pth orbit after given elapsed time
Mp = n*tp; % mean anomaly rad, corresponding to tp th time in ceil(np) orbit
[~,new_True_anomaly] = Eccentric_and_true_anomaly_from_mean_anomaly(Mp,e); % invokes Eccentric_and_true_anomaly_from_mean_anomaly function to calculate
% calculate new eccentric and true anomaly at ceil(np) orbit after elapsed
% time, new_Eccentric_anomaly (rad) , new_True_anomal (rad)
% regression of ascending node
Big_omega_dot = -((3*sqrt(mu_earth)*J2_earth*R_earth^2)/(2*(1-e^2)*a^(7/2)))*cosd(i); % rad/sec regression of ascending node 
Big_omega_dotp = Big_omega + (Big_omega_dot*(180/pi))*delta_t; % deg final regressed RAAN after elapsed time

% advance of perigee
Small_omega_dot = -((3*sqrt(mu_earth)*J2_earth*R_earth^2)/(2*(1-e^2)*a^(7/2)))*(((5/2)*(sind(i))^2)-2); % rad/sec advance of perigee 
Small_omega_dotp = Small_omega + (Small_omega_dot*(180/pi))*delta_t; % deg final advance of perigee after elapsed time

[Perifocal_state_vectors,Geocentric_equatorial_state_vectors] = Geocentric_equi_and_perifocal_state_vector_from_orbital_element(h,e,i,Big_omega_dotp,Small_omega_dotp,(new_True_anomaly)*(180/pi)); % INVOKES Geocentric_equi_and_perifocal_state_vector_from_orbital_element
% function with updates RAAN and argument of perigee to calculate future
% geocentric equatorial and perifocal coordninates.
%% OUTPUT
future_geocentric_equatorial_coordinates = Geocentric_equatorial_state_vectors; % km and km/sec [1X6] future geocentric equatorial coordinates account for earth oblateness
future_perifocal_coordinates = Perifocal_state_vectors; % km and km/sec [1X6] future perifocal coordinates account for earth oblateness
end