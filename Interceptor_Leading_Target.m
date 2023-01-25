function [Delta_v,Time_of_Flight] = Interceptor_Leading_Target(Circular_Orbit_Radius,Angular_Seperation_between_Spacecraft,Target_Revolutions,Interceptor_Revolutions)
% The following code calculate the delta v required and time of flight for
% an interceptor to randevou with the target vehicle.
% Only valid for circular coplanar tranjectories where interceptor is ahead
% of target vehicle and both interceptor and target vehicles are in same
% orbit
% See vallado for further info:-
% REQUIRED INPUT:
% Circular_Orbit_Radius = km Initial circular orbit radius of interceptor and taget
% vehicles
% Angular_Seperation_between_Spacecraft = deg angular seperation between interceptor and target\
% positive if interceptor leads target
% Target_Revolutions = number of revolutions of target vehicle before it Randevous
% Interceptor_Revolutions = number of revolutions of interceptor vehicle before it Randevous
%% Creator:- ANKUR DEVRA 
% Develope Date - 24 March 2022
% Iteration 1 -
%% Starting data
mu_earth = 398600.4418; % Gravitational parameter of earth in km^3/s^2
r = Circular_Orbit_Radius; % km radius of circular orbit
theta = Angular_Seperation_between_Spacecraft*(pi/180); % radian angular seperation between interceptor and target positive if interceptor leads target
k_tgt = Target_Revolutions; % number of revolutions of target vehicle before it Randevous
k_int = Interceptor_Revolutions; % number of revolutions of interceptor vehicle before it Randevous
a_tgt = r; %km semimajor axis of target vehicle
%% Orbit Calculations
Omega_tgt = sqrt(mu_earth/(a_tgt)^3); % rad/sec angular velocity of target vehicle
t_phase = ((k_tgt*(2*pi)) + theta)/Omega_tgt; % sec time of flight for phasing
a_phase = (((mu_earth)*(((t_phase)/(k_int*2*pi))^2)))^(1/3); % km semimajor axis of phasing orbit
%% OUTPUT
Delta_v = 2*(abs(sqrt(((2*mu_earth)/a_tgt) - ((mu_earth)/a_phase)) - sqrt((mu_earth)/(a_tgt)))); %km/sec total delta v required for randevou
Time_of_Flight = t_phase; % sec time of flight till randevou
end