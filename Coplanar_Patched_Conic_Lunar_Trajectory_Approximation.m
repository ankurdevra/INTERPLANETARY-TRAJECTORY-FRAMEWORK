function [Perilune_altitude,TOF_TLI_to_Perilune] = Coplanar_Patched_Conic_Lunar_Trajectory_Approximation(Spacecraft_circular_earth_orbit_altitude,Angular_position_of_spacecraft_relative_to_earth_moon_line_deg,flight_path_angle_at_TLI_deg,Lunar_arrival_angle_deg)
% The following code calculates the perilune altitude and TOF from TLI to
% perilune of a spacecraft in coplanar lunar trajectory using patched
% conics approximation.
% REQUIRED INPUTS:
% Spacecraft_circular_earth_orbit_altitude = km, spacecraft altitude from earth in circular parking orbit
% Angular_position_of_spacecraft_relative_to_earth_moon_line_deg = deg, angular position of spacecraft relative to earth moon line at TLI
% flight_path_angle_at_TLI_deg = deg, flight path angle of spacecraft at TLI
% Lunar_arrival_angle_deg = deg, lunar arrival angle of spacecraft when it enters moon SOI
% OUTPUTS:
% TOF_TLI_to_Perilune = Days, total time of flight from translunar injection to perilune passage of spacecraft
% Perilune_altitude = km, altiude above moon of spacecraft at perilune
% IN OUR ANALYSIS THE ORBIT OF MOON IS ASSUMED TO BE CIRCULAR AROUND EARTH
% OF RADIUS D TO SIMPLIFY OUR ANALYSIS
%% Creator:- ANKUR DEVRA 
% Develope Date - 7 July 2022
% Iteration 1 -
%% Starting data
format long;
mu_earth = 398600.4418; % km^3/s^2 earths gravitaitonal constant
mu_moon = 4904.8695; % km^3/sec^2 moon gravitational parametr
R_earth = 6378.137; % equatorial radius of earth km
R_moon = 1738.1; % equatorial radius of moon km
M_earth = 5.9724*10^(24); % kg, mass of earth
M_moon = 0.07346*10^(24); % kg, mass of moon
D = 384600; % km, moon orbital radius around earth
% IN OUR ANALYSIS THE ORBIT OF MOON IS ASSUMED TO BE CIRCULAR AROUND EARTH
% OF RADIUS D TO SIMPLIFY OUR ANALYSIS
%% INPUTS
r0 = R_earth+Spacecraft_circular_earth_orbit_altitude; % km, spacecraft altitude from earth in circular parking orbit
alpha0 = deg2rad(Angular_position_of_spacecraft_relative_to_earth_moon_line_deg); % rad, angular position of spacecraft relative to earth moon line at TLI
gamma0 = deg2rad(flight_path_angle_at_TLI_deg); % rad, flight path angle of spacecraft at TLI
lambda = deg2rad(Lunar_arrival_angle_deg); % rad, lunar arrival angle of spacecraft when it enters moon SOI
%% CALCULATIONS
R_SOI_moon = D*(M_moon/M_earth)^(2/5); % km, radius of SOI of moon
r0_vec = [-r0*cos(alpha0) -r0*sin(alpha0) 0]; % km, [1X3] position vector of spacecraft at TLI in earth centered non-rotating xyz frame
u_r0 = r0_vec./norm(r0_vec); % [1X3] unit position vector of spacecraft in non-rotating xyz frame with origin of frame centerd at earth
r2_vec = [-R_SOI_moon*cos(lambda) R_SOI_moon*sin(lambda) 0]; % km, [1X3] position vector of spacecraft in non-rotating xyz frame relative to moon (origin of frame at center of moon) when it arrives at SOI of moon
u_r2 = r2_vec./norm(r2_vec); % [1X3] unit position vector of spacecraft in non-rotating xyz frame with origin of frame centered at moon
rm_vec = [D 0 0]; % km, [1X3] vector of position of moon from earth at spacecraft SOI encounter in non-rotating xyz frame with origin of frame centerd at earth
r1_vec = rm_vec + r2_vec; % km, [1X3] position vector of patch point relative to earth of spacecraft in non-rotating xyz frame with origin of frame centerd at earth i.e position of spacecraft wrt to earth at the moment when it arrives at moon SOI
u_r1 = r1_vec./norm(r1_vec); % km, [1X3] unit position vector of patch point relative to earth of spacecraft in non-rotating xyz frame with origin of frame centerd at earth i.e position of spacecraft wrt to earth at the moment when it arrives at moon SOI
delta_theta = acos(dot(u_r0,u_r1)); % rad, sweep angle, difference between the true anomalies of the position vectors r0 and r1
h1 = (sqrt(mu_earth*norm(r0_vec)))*(sqrt((1-cos(delta_theta))/((norm(r0_vec)/norm(r1_vec)) + (sin(delta_theta))*tan(gamma0) - cos(delta_theta)))); % km^2/sec, angular momentum of translunar orbit
% Now calculating lagrange coefficients
f = 1 - (((mu_earth*norm(r1_vec)))/(h1^2))*(1 - cos(delta_theta)); % f lagrange coefficient
g = (((norm(r0_vec))*(norm(r1_vec)))/(h1))*sin(delta_theta); % sec, g lagrange coefficient
g_dot = 1 - (((mu_earth*norm(r0_vec)))/(h1^2))*(1 - cos(delta_theta)); % g dot lagrange coefficient
v0_vec = (1/g).*(r1_vec - f.*r0_vec); % km/sec, [1X3] spacecraft velocity vector at TLI start point
v0_mag = norm(v0_vec); % km/sec, spacecraft velocity at TLI start point
vr_0 = dot(v0_vec,u_r0); % km/sec, radial component of TLI velocity at start to TLI

v1_vec = (1/g).*(g_dot.*r1_vec - r0_vec); % km/sec, [1X3] spacecraft velocity vector when spacecraft reaches moon SOI in TLI trajectory
v1_mag = norm(v1_vec); % km/sec, spacecraft velocity when spacecraft reaches moon SOI in TLI trajectory
vr_1 = dot(v1_vec,u_r1); % km/sec, radial component of spacecraft velocity when it reached moon SOI

e1_vec = (1/mu_earth).*((((v0_mag)^2) - (mu_earth/norm(r0_vec))).*r0_vec - ((norm(r0_vec))*vr_0).*v0_vec); % [1X3] eccentricity vector of Trans lunar tranjectory
e1 = norm(e1_vec); % eccentricity of translunar trajectory

% Now calculating perifocal unit vectors of the translunar trajectory in non-rotating xyz frame with origin of frame centerd at earth
p1_hat = e1_vec./e1; % [1X3] p unit vector of perifocal coords
w1_hat = (cross(r0_vec,v0_vec))./h1; % [1X3] w unit vector of perifocal coords
q1_hat = cross(w1_hat,p1_hat); % [1X3] q unit vector of perifocal coords

a1 = ((h1^2)/(mu_earth))*(1/(1 - e1^2)); % km, semimajor axis of translunar trajectory
T1 = (2*pi)*sqrt((a1^3)/mu_earth); % sec, time period of translunar trajectory

theta0 = acos(dot(p1_hat,u_r0)); % rad, true anomaly of injection point, i.e true anomaly at TLI, which is the angle between perigee of departure trajectory and radial position vector r0 at TLI
% this true anomaly is less that 180 deg, because spacecraft is flying
% outbound towards apogee which lies beyond the patch point

t0 = (T1/(2*pi))*(((2*atan((sqrt((1-e1)/(1+e1)))*tan(theta0/2)))) - e1*sin((2*atan((sqrt((1-e1)/(1+e1)))*tan(theta0/2))))); % sec, time t0 at TLI when spacecraft departs for moon via TLI trajectory

% Patch point is the point when spacecraft enters moon SOI after
% tranvelling via TLI trajectory
theta1 = theta0+delta_theta; % rad, true anomaly of spacecraft at patch point
t1 = (T1/(2*pi))*(((2*atan((sqrt((1-e1)/(1+e1)))*tan(theta1/2)))) - e1*sin((2*atan((sqrt((1-e1)/(1+e1)))*tan(theta1/2))))); % sec, time t1 at patch point when spacecraft arrives at moon via TLI trajectory

delta_t1 = t1-t0; % sec, total time of flight from TLI to SOI (patch point) of moon capture

v_moon = sqrt(mu_earth/D); % km/sec, circular orbital speed of moon
v_moon_vec = [0 v_moon 0]; % km/sec, [1X3] vector of moons velocity at the instant spacecraft crosses moon SOI
omega_moon = v_moon/D; % rad/sec, counterclock wise angular velocity of moon in its assumed circular orbit around earth
% Multiplying this omega_moon with delta_t1 gives moons lead angle. The
% angle throught moon moves as spacecraft flies from TLI to patch point via
% translunar geocentric trajectory

v2_vec = v1_vec - v_moon_vec; % km/sec, [1X3] velocity vector of spacecraft relative to moon
vr_2 = dot(v2_vec,u_r2); % km/sec, radial component of relative velocity of spacecraft at patch point
h2_vec = cross(r2_vec,v2_vec); % km/sec^2 [1X3] angular momentum vector of spacecraft realative to moon at patch point
h2 = norm(h2_vec); % km/sec^2, angular momentum of spacecraft realative to moon at patch point

e2_vec = ((cross(v2_vec,h2_vec))./mu_moon) - u_r2; % [1x3] eccentricity vector of lunar approch trajectory
e2 = norm(e2_vec); % eccentricity of lunar approach trajectory
% if the eccentricity exceeds unity that means inbound orbit is hyperbolic
% realtive to the moon.

p2_hat = (e2_vec)./e2; % [1X3] perifocal unit vector p2 directed from the center of the moon throught the perilune of the approch trajectory
theta2 = (2*pi)-acos(dot(p2_hat,u_r2)); % rad, true anomaly of the patch point on the lunar approach trajectory measured positive clockwise from perilune

t2 = (((h2^3)/((mu_moon^2)*(e2^2 - 1)^(3/2))))*((e2*sinh((2*atanh(((sqrt((e2-1)/(e2+1)))*tan(theta2/2)))))) - (2*atanh(((sqrt((e2-1)/(e2+1)))*tan(theta2/2))))); % sec,time realtive to perilune at patch point.
% according to convention zero time is taken at perigee hence t2 will be
% negative it means it is the time until perilune, starting at sapcecraft arrival at patch point droping towards 0 as spacecraft approach perilune
delta_t2 = 0-t2; % sec, TOF from parch point to perilune

delta_t = delta_t1+delta_t2; % sec, total time of flight from translunar injection to perilune passage of spacecraft
rp_moon = ((h2^2)/(mu_moon))*(1/(1+e2)); % km, perilune radius from center of moon of spacecraft at perilune
vp_moon = (sqrt(1+e2))*(sqrt(mu_moon/rp_moon)); % km/sec, spacecraft speed at perilune realtive to moon
%% OUTPUT
TOF_TLI_to_Perilune = (delta_t)/(3600*24); % days, total time of flight from translunar injection to perilune passage of spacecraft
Perilune_altitude = rp_moon-R_moon; % km, altiude above moon of spacecraft at perilune
end