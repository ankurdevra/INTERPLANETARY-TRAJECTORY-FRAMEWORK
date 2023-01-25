function [Eccentricity,Angular_momentum,Perigee_Velocity,Apogee_Velocity,Semimajor_axis,...
    Period_of_orbit,True_Anomaly_Averaged_Radius] = Orbital_Data_From_Apogee_and_Perigee(Perigee_Altitude,Apogee_Altitude)
% Calculates eccentricity, angular momentum, perigee and apogee velocities,
% semimajor axis, period of orbit and true anomaly averaged radius from
% given input of perigee and apogee altitude for an elliptical orbit around earth. Also plots the variation of
% flight angle with true anomaly.
% REQUIRED INPUTS:
% Perigee_Altitude = Perigee altitude above earth surface in km
% Apogee_Altitude = Apogee altitude above earth surface in km
% OUTPUTS:
% Eccentricity = eccentricity of orbit
% Angular_momentum = angular momentum of orbit km^2/sec
% Perigee_Velocity = velocity at perigee km/sec
% Apogee_Velocity = velocity at apogee km/sec
% Semimajor_axis = semimajor axis of orbit ellips km
% Period_of_orbit = time it takes for one orbit around earth hrs
% True_Anomaly_Averaged_Radius = averages radius of body orbitting earth
%% Creator:- ANKUR DEVRA 
% Develope Date - 23 March 2022
% Iteration 1 -
%% Starting data
format long;
R_earth = 6378.137; % equatorial radius of earth km
mu_earth = 398600.4418; % Gravitational parameter of earth in km^3/s^2
%% Calculations
rp = R_earth + Perigee_Altitude; % radius of perigee in km from center of earth
ra = R_earth + Apogee_Altitude; % radius of apogee in km from center of earth
e = (ra-rp)/(ra+rp); % eccentricity
%h = sqrt(rp*mu_earth*(1-e)) % Angular momentum km^2/sec
% vp = h/rp; % perigee velocity km/sec
% va = h/ra; % apogee velocity km/sec
a = (rp+ra)/2 ;% semimajor axis km
h = sqrt(mu_earth*a*(1-e^2));% Angular momentum km^2/sec
vp = h/rp ;% perigee velocity km/sec
va = h/ra;% apogee velocity km/sec
T = (((2*pi)/sqrt(mu_earth))*a^(3/2))/3600; % time period of orbit hrs
r_theta_bar = sqrt(rp*ra); % true anomaly averaged radius km
%% Outputs
Eccentricity = e;
Angular_momentum = h;
Perigee_Velocity = vp;
Apogee_Velocity = va;
Semimajor_axis = a;
Period_of_orbit = T;
True_Anomaly_Averaged_Radius = r_theta_bar;
%% plotting of Variation of flight path angle with true anomaly for elliptical tranjectory
% all angles in degrees
% figure_number = randi([1 10000],1,2);
% figure(figure_number(1))
% theta = 0:0.01:360;
% gamma = atand((e.*sind(theta))./(1+e.*cosd(theta)));
% plot(theta,gamma,'LineWidth',2);grid on;
% title('Variation of flight path angle with true anomaly in ellipse')
% xlabel('true anomaly (deg)');ylabel('Flight path angle (deg)')
% xlim([0 360]);
end
