function [chase_ellipse_eccentricity,chase_ellipse_angular_momentum,total_delta_v_chase_maneuver] = Chase_maneuver(initial_chaser_and_target_orbit_perigee,initial_chaser_and_target_orbit_apogee,initial_chaser_true_anomaly,initial_target_true_anomaly,desired_time_for_randevou_hrs)
% The following function calculates the total delta-v required to execute a
% chase trajectory as well as eccentricity and angular momentum of chase
% trajectory.
% THE CHASE TRAJECTORY NEED NOT BE A ELLIPSE
% REQUIRED INPUTS:
% initial_chaser_and_target_orbit_perigee = km, perigee of orbit 1 in which both chaser and target lie
% initial_chaser_and_target_orbit_apogee = km, apogee of orbit 1 in which chaser an target lie after chaser randevou with target
% initial_chaser_true_anomaly = deg, initial true anomaly of chaser in orbit 1
% initial_target_true_anomaly = deg, initial true anomaly of target in orbit 1
% desired_time_for_randevou_hrs = hrs, desired time interval at which we want chaser to meet up with traget after exicuting chase manuever
% OUTPUTS:
% chase_ellipse_eccentricity = eccentricity of chase trajectory
% chase_ellipse_angular_momentum = km^2/sec, angular momentum of chase trajectory
% total_delta_v_chase_maneuver = km/sec, total delta-v required to execute the chase manuever
%% Creator:- ANKUR DEVRA 
% Develope Date - 5 July 2022
% Iteration 1 -
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% INPUTS
rp1 = initial_chaser_and_target_orbit_perigee; % km, perigee of orbit 1 in which both chaser and target lie
ra1 = initial_chaser_and_target_orbit_apogee; % km, apogee of orbit 1 in which chaser an target lie after chaser randevou with target
theta_chaser = deg2rad(initial_chaser_true_anomaly); % rad, initial true anomaly of chaser in orbit 1
theta_target = deg2rad(initial_target_true_anomaly); % rad, initial true anomaly of target in orbit 1
delta_t = desired_time_for_randevou_hrs*3600; % sec, desired time interval at which we want chaser to meet up with traget after exicuting chase manuever
%% CALCULATIONS
e1 = (ra1-rp1)/(ra1+rp1); % eccentrity of initial orbit 1 in which both chaser and target lie
h1 = (sqrt(2*mu_earth))*sqrt((ra1*rp1)/(ra1+rp1)); % km^2/sec, angular momentum of initial orbit 1 in which both chaser and target lie

T1 = ((2*pi)/(mu_earth^2))*(h1/(sqrt(1-e1^2)))^3; % sec, time period of orbit 1 in which both chaser and target lie

r_vec_chaser = ((((h1)^2)/(mu_earth))*(1/(1+(e1*cos(theta_chaser))))).*[cos(theta_chaser) sin(theta_chaser) 0]; % km, initial position vector of chaser in perifocal coord in initial orbit 1
v_vec_chaser = (mu_earth/h1).*[-sin(theta_chaser) (e1+cos(theta_chaser)) 0]; % km/sec, initial velocity vector of chaser in perifocal coord in initial orbit 1

E_target = 2*atan(sqrt((1-e1)/(1+e1))*tan(theta_target/2)); % rad, initial eccentric anomaly of target in orbit 1
t_target = (T1/(2*pi))*(E_target - (e1*sin(E_target))); % sec, time since perigee passage of target in orbit 1
t_target_future = t_target + delta_t; % sec, time at which target will be at new position where it will meets with the chaser after elapsed time
M_chaser_future = (2*pi)*(t_target_future/T1); % rad, mean anomaly of target after the elapsed time where it will meet with the chaser

% we invoke the Eccentric_and_true_anomaly_from_mean_anomaly function to
% calculate the eccentric and true anomaly of target after the elapsed time
% where it will meet with the target.
[~,future_true_anomaly_target] = Eccentric_and_true_anomaly_from_mean_anomaly(M_chaser_future,e1);

r_vec_target_future = ((((h1)^2)/(mu_earth))*(1/(1+(e1*cos(future_true_anomaly_target))))).*[cos(future_true_anomaly_target) sin(future_true_anomaly_target) 0]; % km, final position vector of target in perifocal coord in initial orbit 1 where it 
% will meet with the chaser

v_vec_target_future = (mu_earth/h1).*[-sin(future_true_anomaly_target) (e1+cos(future_true_anomaly_target)) 0]; % km/sec, final velocity vector of target in perifocal coord in initial orbit 1 where it 
% will meet with the chaser

% Now we will apply lamberts porblem between initial initial chaser
% position and final chaser position after it has executed chase manuever
% to find the initial velocity of chaser at the start of chase trajectory
% and final velocity of chaser at the end of chase tranjectory
% final chaser position will be the same as final position of target after
% the desired elapsed time

% we invoke Orbital_elements_from_lamberts_problem function to find the
% initial and final state vectors of target while in chase
% trajectory and orbital parameters of chase trajectory
[Angular_momentum,~,Eccentricity,~,~,~,~,~,~,~,~,~,~,~,~,state_vector_at_point_1,state_vector_at_point_2] = Orbital_elements_from_lamberts_problem(r_vec_chaser,r_vec_target_future,desired_time_for_randevou_hrs,'prograde');

v_chaser_in_start_of_chase_trajectory = state_vector_at_point_1(4:6); % km/sec, velocity vector of chaser at the start of chase trajectory
v_chaser_in_end_of_chase_trajectory = state_vector_at_point_2(4:6); % km/sec, velocity vector of chaser at the end of chase trajectory

delta_v_chaser_start = norm(v_chaser_in_start_of_chase_trajectory-v_vec_chaser); % km/sec, delta-v required by chaser to initiate a chase trajectory
delta_v_chaser_end = norm(v_vec_target_future-v_chaser_in_end_of_chase_trajectory); % km/sec, delta-v required by chaser to end the chase trajectory
%% OUTPUT
chase_ellipse_eccentricity = Eccentricity; % eccentricity of chase trajectory
chase_ellipse_angular_momentum = Angular_momentum; % km^2/sec, angular momentum of chase trajectory
total_delta_v_chase_maneuver = delta_v_chaser_start+delta_v_chaser_end; % km/sec, total delta-v required to execute the chase manuever
end