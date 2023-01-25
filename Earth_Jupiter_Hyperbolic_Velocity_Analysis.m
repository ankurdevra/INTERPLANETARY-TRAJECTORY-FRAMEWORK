%% Jupiter flyby analysis
function [departure_hyperbolic_excess_velocity,arrival_hyperbolic_excess_velocity,departure_hyperbolic_velocity_vec,arrival_hyperbolic_velocity_vec] = Earth_Jupiter_Hyperbolic_Velocity_Analysis(SOI_Earth_Departure_Date_and_Time,SOI_Jupiter_Arrival_Date_and_Time)
format long
mu_jupiter = 126686531.900;% km^3/sec^2 gravitational prameter of jupiter
R_jupiter_radius = 71492; % km Jupiter Equatorial Radius
Jupiter_state_vectors = readmatrix("Jupiter_state_vectors_1_jan_2000_to_31_dec_2099.txt");
Julian_dates = rmmissing(Jupiter_state_vectors(:,1)); % julian dates from 1 sept 2036 to 31 dec 2039
% JUPITER STATE VECTOR DATA
X_Jupiter = rmmissing(Jupiter_state_vectors(:,3)); % km, x position of jupiter vector
Y_Jupiter = rmmissing(Jupiter_state_vectors(:,4)); % km, y position of jupiter vector
Z_Jupiter = rmmissing(Jupiter_state_vectors(:,5)); % km, z position of jupiter vector
VX_Jupiter = rmmissing(Jupiter_state_vectors(:,6)); % km, vx velocity of jupiter vector
VY_Jupiter = rmmissing(Jupiter_state_vectors(:,7)); % km, vy velocity of jupiter vector
VZ_Jupiter = rmmissing(Jupiter_state_vectors(:,8)); % km, vz velocity of jupiter vector
% EARTH STATE VECTOR DATA
Earth_state_vectors = readmatrix("Earth_state_vectors_1_jan_2000_to_31_dec_2099.txt");
X_Earth = rmmissing(Earth_state_vectors(:,3)); % km, x position of Earth vector
Y_Earth = rmmissing(Earth_state_vectors(:,4)); % km, y position of Earth vector
Z_Earth = rmmissing(Earth_state_vectors(:,5)); % km, z position of Earth vector
VX_Earth = rmmissing(Earth_state_vectors(:,6)); % km, vx velocity of Earth vector
VY_Earth = rmmissing(Earth_state_vectors(:,7)); % km, vy velocity of Earth vector
VZ_Earth = rmmissing(Earth_state_vectors(:,8)); % km, vz velocity of Earth vector
%
Departure_date = SOI_Earth_Departure_Date_and_Time;%[2036 9 1 0 0 0]; % [yyyy mm dd] departure date and UT time from earth, when spacecraft leave earth SOI
Departure_date_julian = juliandate(Departure_date); % departure date from earth julian date
Arrival_date = SOI_Jupiter_Arrival_Date_and_Time;%[2039 1 1 0 0 0]; % [yyyy mm dd] arrival date and UT time to jupiter, when soacecraft reached jupiter SOI
Arrival_date_julian = juliandate(Arrival_date); %  arrival date to jupiter julian date
[row_departure_earth,column_departure_earth] = find(Julian_dates==Departure_date_julian); % julian date index of departure
[row_arrival_jupiter,column_arrival_jupiter] = find(Julian_dates==Arrival_date_julian); % julian date index of arrival

Earth_departure_state_vector = [X_Earth(row_departure_earth,column_departure_earth) Y_Earth(row_departure_earth,column_departure_earth) Z_Earth(row_departure_earth,column_departure_earth) VX_Earth(row_departure_earth,column_departure_earth) VY_Earth(row_departure_earth,column_departure_earth) VZ_Earth(row_departure_earth,column_departure_earth)]; % km, km/s state vector of earth at the time of departure
Jupiter_arrival_state_vector = [X_Jupiter(row_arrival_jupiter,column_arrival_jupiter) Y_Jupiter(row_arrival_jupiter,column_arrival_jupiter) Z_Jupiter(row_arrival_jupiter,column_arrival_jupiter) VX_Jupiter(row_arrival_jupiter,column_arrival_jupiter) VY_Jupiter(row_arrival_jupiter,column_arrival_jupiter) VZ_Jupiter(row_arrival_jupiter,column_arrival_jupiter)]; % km, km/s state vector of jupiter at the time of arrival
TOF_hrs = (Arrival_date_julian-Departure_date_julian)*24; % time of flight from earth to jupiter in hrs
norm_earth_posistion = norm(Earth_departure_state_vector(1:3));
norm_earth_velocity = norm(Earth_departure_state_vector(4:6));
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,state_vector_at_point_1,state_vector_at_point_2] = Heliocentric_Orbital_elements_from_lamberts_problem(Earth_departure_state_vector(1:3),Jupiter_arrival_state_vector(1:3),TOF_hrs,'prograde');

departure_hyperbolic_velocity_vec = state_vector_at_point_1(4:6) - Earth_departure_state_vector(4:6); % km/sec hyperbolic departure excess velocity vector of spacecraft in ICRF coord frame
arrival_hyperbolic_velocity_vec = state_vector_at_point_2(4:6) - Jupiter_arrival_state_vector(4:6); % km/sec hyperbolic arrival excess velocity vector of spacecraft in ICRF coord frame
departure_hyperbolic_excess_velocity = norm(departure_hyperbolic_velocity_vec); % km/sec hyperbolic departure excess velocity of spacecraft to complete travel in desired given amount of time
arrival_hyperbolic_excess_velocity = norm(arrival_hyperbolic_velocity_vec); % km/sec hyperbolic arrival excess velocity of spacecraft to complete travel in desired given amount of time

end
%[departure_hyperbolic_excess_velocity_,arrival_hyperbolic_excess_velocity_,TOF_hrs_]  = Interplanetary_Trajectory('Earth','Jupiter',[2036 9 1 0 0 0],[2039 1 1 0 0 0]);



