function [total_delta_v_for_Apse_line_rotation,true_anomaly_required_for_Apse_line_rotation,orientation_of_delta_v_required_for_Apse_line_rotation] = Apse_line_rotation_from_orbit1_to_orbit2(initial_orbit_perigee,initial_orbit_apogee,final_orbit_perigee,final_orbit_apogee,desired_apse_line_rotation)
% The following function calculates the total delta-v required, true
% anomaly required and orientation of total delta-v velocity vector to
% achieve the desired apse line rotation from initial orbit 1 to final
% orbit 2
% REQUIRED INPUTS:
% initial_orbit_perigee = km, perigee of initial orbit 1
% initial_orbit_apogee = km, apogee of initial orbit 1
% final_orbit_perigee = km, perigee of final desired orbit 2
% final_orbit_apogee = km, apogee of final desired orbit 2
% desired_apse_line_rotation = deg, desired apse line rotation to go from initial to final orbit
% OUTPUTS:
% total_delta_v_for_Apse_line_rotation = km/sec, delta-v required to 
% achieve the desired apse line rotation from initial orbit 1 to final
% orbit 2
% true_anomaly_required_for_Apse_line_rotation = deg, true anomaly of body in orbit 1 to achieve the desire
% apse line rotation from initial orbit 1 to final
% orbit 2
% orientation_of_delta_v_required_for_Apse_line_rotation = deg, angle phi vector of delta-v need to make with 
% local horizon to achieve the desired apse line rotation from initial orbit 1 to final
% orbit 2
%% Creator:- ANKUR DEVRA 
% Develope Date - 5 July 2022
% Iteration 1 -
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% INPUTS
rp1 = initial_orbit_perigee; % km, perigee of initial orbit 1
ra1 = initial_orbit_apogee; % km, apogee of initial orbit 1
rp2 = final_orbit_perigee; % km, perigee of final desired orbit 2
ra2 = final_orbit_apogee; % km, apogee of final desired orbit 2
eta = deg2rad(desired_apse_line_rotation); % radian, apse line rotation to go from initial to final orbit 
%% CALCULATIONS
e1 = (ra1-rp1)/(ra1+rp1); % eccentricity of initial orbit
e2 = (ra2-rp2)/(ra2+rp2); % eccentricity of final orbit
h1 = (sqrt(2*mu_earth))*sqrt((rp1*ra1)/(rp1+ra1)); % km^2/sec angular momentum of initial orbit 1
h2 = (sqrt(2*mu_earth))*sqrt((rp2*ra2)/(rp2+ra2)); % km^2/sec angular momentum of final orbit 2
% intermediate values for calculations
a = (e1*(h2)^2) - ((e2*(h1)^2)*cos(eta)); % km^4/sec^2
b = -((e2*(h1)^2)*sin(eta)); % km^4/sec^2
c = (h1)^2 - (h2)^2; % km^4/sec^2
phi = atan(b/a); % rad
% end of intermediate calculation

% The initial orbit intersects the final orbit at two points I and J.
% we will use the point with smaller true anomaly for apse line rotation
% manuver
theta1 = (phi + acos((c/a)*cos(phi))); % rad, true anomaly of intersection point 1 (I)
theta2 = (phi - acos((c/a)*cos(phi))); % rad, true anomaly of intersection point 2 (J)
if theta2<0
    theta2 = 2*pi + theta2; % rad, makes sure that true anomaly is positive quantity
end
% compare which is small 
if theta1<theta2
    theta_1 = theta1; % rad, true anomaly of manuvering point (I) from initial orbit 1
else
    theta_1 = theta2; % rad, true anomaly of manuvering point (I) from initial orbit 1
end
theta_2 = theta_1-eta; % rad, true anomaly of manuvering point (I) from final orbit 2
r = (((h1)^2)/(mu_earth))*(1/(1+e1*cos(theta_1))); % km, radian coordinate of manuvering point I
% calculating velocity component and flight path angle for initial orbit 1
% at point I
v_t1 = h1/r; % km/sec, transverse component of velocity at point I in initial orbit 1
v_r1 = (mu_earth/h1)*e1*sin(theta_1); % km/sec, radial component of velocity at point I in initial orbit 1
gamma1 = atan(v_r1/v_t1); % rad, flight path angle at point I in intial orbit 1
v1 = sqrt((v_r1)^2 + (v_t1)^2); % km/sec, speed of body at point I in initial orbit 1

% calculating velocity component and flight path angle for final orbit 2
% at point I
v_t2 = h2/r; % km/sec, transverse component of velocity at point I in final orbit 2
v_r2 = (mu_earth/h2)*e2*sin(theta_2); % km/sec, radial component of velocity at point I in final orbit 2
gamma2 = atan(v_r2/v_t2); % rad, flight path angle at point I in final orbit 2
v2 = sqrt((v_r2)^2 + (v_t2)^2); % km/sec, speed of body at point I in final orbit 2

phi = atan((v_r2-v_r1)/(v_t2-v_t1)); % rad, angle phi vector of delta-v need to make with 
% local horizon to achieve the desired apse line rotation from initial orbit 1 to final
if phi<0
    phi = pi + phi; % makes sure phi is positive 
end
%% OUTPUT
total_delta_v_for_Apse_line_rotation = sqrt((v1)^2 + (v2)^2 - (2*v1*v2)*cos(gamma2-gamma1)); % km/sec, delta-v required to 
% achieve the desired apse line rotation from initial orbit 1 to final
% orbit 2
true_anomaly_required_for_Apse_line_rotation = rad2deg(theta_1); % deg, true anomaly of body in orbit 1 to achieve the desire
% apse line rotation from initial orbit 1 to final
% orbit 2
orientation_of_delta_v_required_for_Apse_line_rotation = rad2deg(phi); % deg, angle phi vector of delta-v need to make with 
% local horizon to achieve the desired apse line rotation from initial orbit 1 to final
% orbit 2
end