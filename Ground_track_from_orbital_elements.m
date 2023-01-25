function [Ground_track_plot] = Ground_track_from_orbital_elements(Angular_momentum,eccentricity,inclination,RAAN,Argument_of_perigee,True_anomaly,elapsed_time_hrs)
% The following function plots the ground track of a sattelite orbiting
% earth from its orbital elements till the specified elpased time. It also
% accounts for the oblateness of earth and rotation of earth to compute the
% ground track.
% REQUIRED INPUTS:
% Angular_momentum = angular momentum of orbiting body around earth km^2/sec
% eccentricity = orbital eccentricity
% inclination  = orbital inclination deg
% RAAN = Right ascension of ascending node deg
% Argument_of_perigee = argument of perigee deg
% True_anomaly = true anomaly deg
% elapsed_time_hrs = elapsed time till which we want ground track.
% OUTPUT
% Ground_track_plot = GROUND track of sattelite orbiting earth
%% Creator:- ANKUR DEVRA 
% Develope Date - 3 July 2022
% Iteration 1 -
%% CALCULATIONS
time = 0:0.001:elapsed_time_hrs;  % time vector
future_right_ascension = zeros(1,length(time)); % initializing matrix
future_decilantion = zeros(1,length(time)); % initializing matrix

% calls future_RA_Dec_relative_rotating_earth_from_orbital_elements
% function at each time step and calculate RA and Dec.
for i = 1:length(time)
    [future_right_ascension(i+1),future_decilantion(i+1)] = future_RA_Dec_relative_rotating_earth_from_orbital_elements(Angular_momentum,eccentricity,inclination,RAAN,Argument_of_perigee,True_anomaly,time(i));
end
future_right_ascension_vector = future_right_ascension; % stores the values of RA
future_decilantion_vector = future_decilantion; % stores the values of Dec 
% Plots RA and Dec values on the surface of earth.
figure(1)
hold on;
earth_pic = imread('2_no_clouds_4k.jpg');
Fit = earth_pic(end:-1:1,:,:);
image([0 360],[-90 90],Fit)
Ground_track_plot = plot(future_right_ascension_vector,future_decilantion_vector,'k.');
axis([0,360,-90,90])
daspect([1 1 1]);
xlabel('Right Ascension (deg)');ylabel('Declination (deg)')
title('Ground Track of Satellite')
end