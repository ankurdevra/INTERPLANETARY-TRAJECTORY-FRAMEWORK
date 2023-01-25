function [minimum_delta_v,periapsis_radius,aiming_radius,angle_bewtween_periapsis_and_mars_velocity_vector] = hohmann_transer_from_earth_to_mars_capture_data(time_period_of_desired_orbit_around_mars_hrs)
% The following function calculates the minimum delta-v required to place
% a spacecraft in the desired orbit with desired time period around mars
% after hohmann interplanetary journey from earth to mars along with periapses
% radius, aiming radius and angle between periapses and mars velocity
% vector for the spacecraft to be places in desired orbit
% INPUTS:
% time_period_of_desired_orbit_around_mars_hrs = sec, desired time period of orbit around mars
% OUTPUTS:
% minimum_delta_v = km/sec, minimum delta-v required to place the spacecraft in desired orbit around mars
% periapsis_radius = km, periapses radius of orbit around mars for the given hyperbolic excess speed
% aiming_radius = km, aiming radius of spacecraft to get into desired orbit around mars
% angle_bewtween_periapsis_and_mars_velocity_vector = deg, location of periapses wrt apse line and arrival hyperbola
% SIMPLIFIED ANALYSIS ASSUMING ORBITS OF MARS AND EARTH ARE CIRCULAR AND
% COPLANAR
%% Creator:- ANKUR DEVRA 
% Develope Date - 7 July 2022
% Iteration 1 -
%% Starting data
R_earth = 149.6*10^(6); % km Earth orbital Radius
R_mars_radius = 3396.2; % km mars Equatorial Radius
R_mars = 227.9*10^(6); % km Mars orbital radius
mu_mars = 42828.37; %km^3/s^2 mars gravitaitonal parameter
mu_sun = 1.32712440018*10^(11); % km^3/s^2 gravitational parameter
%% INPUTS
T = 3600*time_period_of_desired_orbit_around_mars_hrs; % sec, desired time period of orbit around mars
%% CALCULATIONS
v_inf = sqrt(mu_sun/R_mars)*(1-sqrt((2*R_earth)/(R_mars+R_earth))); % km/sec, hyperbolic excess speed spacecraft must arrive at mars to be captured in its SOI
a = ((T*sqrt(mu_mars))/(2*pi))^(2/3); % km, semimajor axis of the orbit around mars for the given desired orbital time period of spacecraft around mars
e = (2*mu_mars)/(a*v_inf^2)-1; % eccentricity of the orbit around mars for the given desired orbital time period of spacecraft around mars
rp = ((2*mu_mars)/(v_inf)^2)*((1-e)/(1+e)); % km, periapses radius of orbit around mars for the given hyperbolic excess speed
beta = acos(1/(1+((rp*v_inf^2)/mu_mars))); % rad, Beta gives the orientation
% of the apse line of the arrival hyperbola to the planet's heliocentric
% velocity vector
%% OUTPUT
minimum_delta_v = v_inf*sqrt((1-e)/2); % km/sec, minimum delta-v required to place the spacecraft in desired orbit around mars
periapsis_radius = rp; % km, periapses radius of orbit around mars for the given hyperbolic excess speed
aiming_radius = rp*sqrt(2/(1-e)); % km, aiming radius of spacecraft to get into desired orbit around mars
angle_bewtween_periapsis_and_mars_velocity_vector = rad2deg(beta); % deg, location of periapses wrt apse line and arrival hyperbola
end