function [delta_v_earth_mars_trajectory,location_of_perigee,amount_of_propellant_as_percent_of_spacecraft] = earth_to_mars_interplanetary_trajectory_delta_v_and_perigee(radius_of_circular_parking_orbit,specific_impulse_of_spacecraft)
% The following code claculates the delta-v, location of perigee (wrt apse
% line and planet's heliocentric velocity vector) and propellant as a
% percentage of spacecraft mass required for hyperbolic departure from
% earth to mars using hohmann transfer.
% REQUIRED INPUTS:
% radius_of_circular_parking_orbit = km, radius of initial circular parking orbit (above earth surface), from which spacecraft execute interplanetary manuever
% specific_impulse_of_spacecraft = sec, specific impulse of spacecraft engine
% OUTPUTS:
% delta_v_earth_mars_trajectory = km/sec delta-v required to put the spacecraft onto the hyperbolic departure trajectory
% location_of_perigee = deg,location of periapsis where delta_v maneuver must take place to put the 
% spacecraft in hyperbolic departure trajectory.
% amount_of_propellant_as_percent_of_spacecraft = percent (%), prior to delta-v
% maneuver percentage of spacecraft mass that must be propellant.
% SIMPLIFIED ANALYSIS ASSUMING ORBITS OF MARS AND EARTH ARE CIRCULAR AND
% COPLANAR
%% Creator:- ANKUR DEVRA 
% Develope Date - 7 July 2022
% Iteration 1 -
%% Starting data
R_earth = 149.6*10^(6); % km Earth orbital Radius
R_earth_radius = 6378.137; % km Earth Equatorial Radius
R_mars = 227.9*10^(6); % km Mars orbital radius
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal parameter
mu_sun = 1.32712440018*10^(11); % km^3/s^2 sun gravitational parameter
g0 = 9.81/1000; % km/sec^2, sea level gravitational acceleration of earth
%% INPUTS
rp = radius_of_circular_parking_orbit; % km, radius of initial circular parking orbit, from which spacecraft execute interplanetary manuever
Isp = specific_impulse_of_spacecraft; % sec, specific impulse of spacecraft engine
%% CALCULATIONS
v_inf = sqrt(mu_sun/R_earth)*(sqrt(((2*R_mars)/(R_earth+R_mars))) - 1); % km/sec, hyperbolic excess speed relative to home planet
% of the departure hyperbola. This
% is the velocity at which the spacecraft must arrive at the sphere of
% influence of the home planet to escape the home planet gravitational pull
% and travel to the target planet
v_c = sqrt((mu_earth)/(R_earth_radius+rp)); % km/sec, speed of spacecraft in circular orbit around earth
v_p = sqrt((v_inf)^2 + (2*mu_earth)/(R_earth_radius+rp)); % km/sec, speed required at the periapsis of circular parking orbit to exit the
% sphere of influence of earth with the required hyperbolic excess speed
e = 1+((((rp+R_earth_radius)*(v_inf)^2))/mu_earth); % eccentricity of departure hyperbola
beta = acos(1/e); % rad, location of periapsis where delta_v maneuver must take place to put the 
% spacecraft in hyperbolic departure trajectory. Beta gives the orientation
% of the apse line of the departure hyperbola to the planet's heliocentric
% velocity vector
%% OUTPUT
delta_v_earth_mars_trajectory = v_p-v_c; % km/sec delta-v required to put the spacecraft onto the hyperbolic departure trajectory
location_of_perigee = rad2deg(beta); % deg, location of periapsis where delta_v maneuver must take place to put the 
% spacecraft in hyperbolic departure trajectory.
amount_of_propellant_as_percent_of_spacecraft = (1-exp(-delta_v_earth_mars_trajectory/(Isp*g0)))*100; % percent (%), prior to delta-v
% maneuver percentage of spacecraft mass that must be propellant.
end