function [departure_hyperbolic_excess_velocity,arrival_hyperbolic_excess_velocity,TOF_hrs]  = Interplanetary_Trajectory(Departure_planet,Arrival_planet,SOI_Departure_date_and_time,SOI_Arrival_date_and_time)
% The following function calculates the required departure and arrival hyperbolic excess velocity to travel
% from departure planet to arrival planet in given specfied amount of time.
% REQUIRED INPUTS:
% Departure_planet = departure planet name 'Earth','Mars' etc
% Arrival_planet = arrival planet name 'Mars', 'Saturn' etc
% SOI_Departure_date_and_time = [YYYY MM D UT_HRS UT_MIN UT_SEC],
% year,month,date along with UT hrs,min,sec when the spacecraft exit the
% SOI of departure planet
% SOI_Arrival_date_and_time = [YYYY MM D UT_HRS UT_MIN UT_SEC],
% year,month,date along with UT hrs,min,sec when the spacecraft enters the
% SOI of arrival planet
% OUTPUTS:
% departure_hyperbolic_excess_velocity = km/sec hyperbolic departure excess
% velocity of spacecraft to complete travel in desired given amount of time
% arrival_hyperbolic_excess_velocity = km/sec hyperbolic arrival excess velocity of spacecraft 
% to complete travel in desired given amount of time
% *ALL VECTORS ARE WRT ICRF COORD FRAME*
%% Creator:- ANKUR DEVRA 
% Develope Date - 8 July 2022
% Iteration 1 -
%% CALCULATIONS
JD_departure = juliandate(SOI_Departure_date_and_time(1),SOI_Departure_date_and_time(2),SOI_Departure_date_and_time(3),SOI_Departure_date_and_time(4),SOI_Departure_date_and_time(5),SOI_Departure_date_and_time(6)); % Julian day number from 
% date and UT time of spacecraft when it exits the departure planet SOI.
JD_arrival = juliandate(SOI_Arrival_date_and_time(1),SOI_Arrival_date_and_time(2),SOI_Arrival_date_and_time(3),SOI_Arrival_date_and_time(4),SOI_Arrival_date_and_time(5),SOI_Arrival_date_and_time(6));% Julian day number from 
% date and UT time of spacecraft when it enters the arrival planet SOI.
[R_departure_planet,V_departure_planet] = planetEphemeris(JD_departure,'Sun',Departure_planet,'432t'); % km and km/sec,heliocentric coordinates of departure planet at the
% time of departure in ICRF coords frame
[R_arrival_planet,V_arrival_planet] = planetEphemeris(JD_arrival,'Sun',Arrival_planet,'432t');% km and km/sec,heliocentric coordinates of arrival planet at the
% time of arrival in ICRF coords frame
R_departure_spacecraft = R_departure_planet; % km, heliocentric position of spacecraft at the time of SOI departure from departure planet in ICRF coords frame
R_arrival_spacecraft = R_arrival_planet; % km, heliocentric position of spacecraft at the time of SOI arrival of arrival planet in ICRF coords frame
TOF_hrs = (JD_arrival-JD_departure)*24; % hrs, time of flight in hrs from departure planet to arrival planet
% Now solving lamberts problem between SPACECRAFT departure position and
% spacecraft arrival position
% using the TOF to find the required departure velocity and arrival velocity of
% spacecraft between the desired TOF
% *THE OUTPUT STATE VECTORS AND COORDINATE DEPENDENT ORBITAL DATA WILL BE WITH RESPECT TO ICRF
% FRAME*
% eccentricity,angular momentum,semimajor axis,semilatus rectum etc of an orbit are independent of coordinate
% system used as they are physical properties of an orbit and does not
% require a reference point
% ICRF COORDS AND HELIOCENTRIC COORDS DO NOT ALIGHT WITH EACH OTHER
% ORIGIN OF ICRF SYSTEM IS AT BARYCENTER OF OUR SOLARSYSTEM
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,state_vector_at_point_1,state_vector_at_point_2] = Heliocentric_Orbital_elements_from_lamberts_problem(R_departure_spacecraft,R_arrival_spacecraft,TOF_hrs,'prograde');
departure_hyperbolic_velocity_vec = state_vector_at_point_1(4:6) - V_departure_planet; % km/sec hyperbolic departure excess velocity vector of spacecraft in ICRF coord frame
arrival_hyperbolic_velocity_vec = state_vector_at_point_2(4:6) - V_arrival_planet; % km/sec hyperbolic arrival excess velocity vector of spacecraft in ICRF coord frame
%% OUTPUT
departure_hyperbolic_excess_velocity = norm(departure_hyperbolic_velocity_vec); % km/sec hyperbolic departure excess velocity of spacecraft to complete travel in desired given amount of time
arrival_hyperbolic_excess_velocity = norm(arrival_hyperbolic_velocity_vec); % km/sec hyperbolic arrival excess velocity of spacecraft to complete travel in desired given amount of time
end