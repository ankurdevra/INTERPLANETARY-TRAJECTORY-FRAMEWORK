function [chaser_relative_position_wrt_target,chaser_relative_velocity_wrt_target,chaser_relative_accelaration_wrt_target] = chaser_relative_to_target_vectors_from_orbital_elements(chaser_orbit_angular_momentum,chaser_orbit_eccentricity,chaser_orbit_inclination,chaser_orbit_RAAN,chaser_orbit_argument_of_perigee,chaser_orbit_true_anomaly,target_orbit_angular_momentum,target_orbit_eccentricity,target_orbit_inclination,target_orbit_RAAN,target_orbit_argument_of_perigee,target_orbit_true_anomaly)
% The following function calculates the relative position, velocity and
% acceleration of chaser relative to target from their orbital elements. That is motion of chaser as
% seen from the targets perseptive.
% REQUIRED INPUTS:
% chaser_orbit_angular_momentum = km/sec^2, angular momentum of chaser orbit
% chaser_orbit_eccentricity = eccentricity of chaser orbit
% chaser_orbit_inclination = deg, orbital inclination of chaser orbit
% chaser_orbit_RAAN = deg, RAAN of chaser orbit
% chaser_orbit_argument_of_perigee = deg, argument of perigee of chaser orbit
% chaser_orbit_true_anomaly = deg, true anomaly of chaser spacecraft
% target_orbit_angular_momentum = km/sec^2, angular momentum of target orbit
% target_orbit_eccentricity = eccentricity of target orbit
% target_orbit_inclination = deg, orbital inclination of target orbit
% target_orbit_RAAN = deg, RAAN of target orbit
% target_orbit_argument_of_perigee = deg, argument of perigee of target orbit
% target_orbit_true_anomaly = deg, true anomaly of target spacecraft
% OUTPUTS:
% chaser_relative_position_wrt_target =  [1X3] km, relative position of chaser as seen from target
% chaser_relative_velocity_wrt_target =  [1X3] km/sec relative velocity of chaser as seen from target
% chaser_relative_accelaration_wrt_target = [1X3] km/sec^2 relative acceleration of chaser as seen from target
%% Creator:- ANKUR DEVRA 
% Develope Date - 6 July 2022
% Iteration 1 -
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% INPUTS
% invokes Geocentric_equi_and_perifocal_state_vector_from_orbital_element
% function to calculate the current geocentric state vectors of chaser and
% target spacecraft in earths orbit.
% GEOCENTRIC EQUATORIAL FRAME IS INERTIAL FRAME
% GEOCENTRIC [I J K]; LVLH [i j k]
[~,Chaser_Geocentric_equatorial_state_vectors] = Geocentric_equi_and_perifocal_state_vector_from_orbital_element(chaser_orbit_angular_momentum,chaser_orbit_eccentricity,chaser_orbit_inclination,chaser_orbit_RAAN,chaser_orbit_argument_of_perigee,chaser_orbit_true_anomaly); % km and km/sec, current geocentric state vectors of chaser spacecraft
[~,Target_Geocentric_equatorial_state_vectors] = Geocentric_equi_and_perifocal_state_vector_from_orbital_element(target_orbit_angular_momentum,target_orbit_eccentricity,target_orbit_inclination,target_orbit_RAAN,target_orbit_argument_of_perigee,target_orbit_true_anomaly); % km and km/sec, current geocentric state vectors of target spacecraft
r_chaser_vec = Chaser_Geocentric_equatorial_state_vectors(1:3); % km, [1X3] current geocentric position vector of chaser spacecraft
r_chaser_mag = norm(r_chaser_vec); % km, current geocentric position of chaser spacecraft
v_chaser_vec = Chaser_Geocentric_equatorial_state_vectors(4:6); % km/sec, [1X3] current geocentric velocity vector of chaser spacecraft
v_chaser_mag = norm(v_chaser_vec); % km/sec, current geocentric velocity of chaser spacecraft
r_target_vec = Target_Geocentric_equatorial_state_vectors(1:3); % km, [1X3] current geocentric position vector of target spacecraft
r_target_mag = norm(r_target_vec); % km, current geocentric position of target spacecraft
v_target_vec = Target_Geocentric_equatorial_state_vectors(4:6); % km/sec, [1X3] current geocentric velocity vector of target spacecraft
v_target_mag = norm(v_target_vec); % km/sec, current geocentric velocity of target spacecraft
%% CALCULATIONS
h_target_vec = cross(r_target_vec,v_target_vec); % km^2/sec, [1X3] current geocentric angular momentum vector of target spacecraft
h_target_mag = norm(h_target_vec); % km^2/sec, current geocentric angular momentum of target spacecraft
% Now calculating the unit vectors i,j,k of the comoving LVLH frame
% LVLH FRAME IS NON INERTIAL
i_hat = r_target_vec./r_target_mag; % [1X3] i unit vector of LVLH frame in terms of geocentric coords
k_hat = h_target_vec./h_target_mag; % [1X3] k unit vector of LVLH frame in terms of geocentric coords
j_hat = cross(k_hat,i_hat); % [1X3] j unit vector of LVLH frame in terms of geocentric coords
% We now form the [Q]_X_x orthogonal DCM. This is used to convert inertial
% frame of reference to non inertial frame of reference.
Q_X_x = [i_hat(1:3);j_hat(1:3);k_hat(1:3)]; % [3X3] DCM to go from inertial geocentric equatorial frame to 
% non inertal LVLH frame. Row 1 are components of i unit vectors; Row 2 are
% components of j unit vectors; Row 3 are components of k unit vectors.
omega = h_target_vec./((r_target_mag)^2); % [1X3] rad/sec, angular velocity vector of LVLH xyz axes attached to target in geocentric coords
omega_dot = -2.*(((dot(v_target_vec,r_target_vec))./((r_target_mag)^2)).*omega); % [1X3] rad/sec^2 angular accelaration vector  of LVLH xyz axes in geocentric coords
absolute_acceleration_target_vec = ((-mu_earth)/((r_target_mag)^3)).*(r_target_vec); % [1X3] km/sec^2 absolute acceleration of target spacecraft in geocentric coords
absolute_acceleration_chaser_vec = ((-mu_earth)/((r_chaser_mag)^3)).*(r_chaser_vec); % [1X3] km/sec^2 absolute acceleration of target spacecraft in geocentric coords

relative_position_vec = r_chaser_vec-r_target_vec; % [1X3] km, chaser relative to target position vector in geocentric coords
relative_velocity_vec = v_chaser_vec-v_target_vec-cross(omega,relative_position_vec); % [1X3] km/sec chaser relative to target velocity vector in geocentric coords
relative_acceleration_vec = absolute_acceleration_chaser_vec-absolute_acceleration_target_vec-cross(omega_dot,relative_position_vec)-cross(omega,cross(omega,relative_position_vec))-2.*(cross(omega,relative_velocity_vec)); % [1X3] km/sec^2 chaser relative to target acceleration vector in geocentric coords
% The above vector calculation was done wrt geocentric equatorial coords
% hence all of the above vectors have geocentric basis vectors namely [I J K]
% now we will use the DCM to convert geocentric to LVLH [I J K] ---> [i j k]
% we will go from inertial geocentric to non inertial LVLH frame
%% OUTPUT
chaser_relative_position_wrt_target = (Q_X_x*relative_position_vec')'; % [1X3] km, relative position of chaser as seen from target
chaser_relative_velocity_wrt_target = (Q_X_x*relative_velocity_vec')'; % [1X3] km/sec relative velocity of chaser as seen from target
chaser_relative_accelaration_wrt_target = (Q_X_x*relative_acceleration_vec')'; % [1X3] km/sec^2 relative acceleration of chaser as seen from target
%% PLOTTING THE MOTION OF CHASER RELATIVE TO TARGET
% INCOMPLETE FOR LOOP, NEEDS TO BE COMPLETED LATER IN SOME OTHER m file

% T_target = ((2*pi)/((mu_earth)^2))*((h_target_mag/(sqrt(1-target_orbit_eccentricity^2))))^3; % sec, time period of target spacecraft orbit
% t0 = 0; % sec, initial starting time for plot
% m = 1; % no of target orbit to plot the motion, currently set 1. means how will the motion look for one orbit of target
% tf = t0 + m*T_target; % sec, final ending time for plot
% steps = 100; % amount of steps or points to be taken for plotting
% dt = tf/steps; % sec, amount to move forward in time for plotting
% for i = 1:tf
%     t=t+dt; % sec, increment of time per iteration
%     [target_state_vector,~] = future_state_vector_from_initial_using_universal_anomaly(Target_Geocentric_equatorial_state_vectors,t); % km and km/sec [1X6] target geocentric state vector after elapsed time t
%     [chaser_state_vector,~] = future_state_vector_from_initial_using_universal_anomaly(Chaser_Geocentric_equatorial_state_vectors,t); % km and km/sec [1X6] chaser geocentric state vector after elapsed time t
% end


end