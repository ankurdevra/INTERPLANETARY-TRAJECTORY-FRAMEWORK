function [State_vectors_JPL_Horizons] = State_Vector_of_Planets_at_specified_Date_and_Time(State_Vector_date_and_time,Planet)
% The following code extracts the state vectors of a specified planet.
% The state vectors are obtained from NASA JPL HORIZONS with reference
% frame Ecliptic of J2000.0
% *NOTE - CURRENTLY ONLY EXTRACTS STATE VECTORS AT SPECIFIED JULIAN DATE
% ONLY, ONLY JULIAN DATES ENDING IN XXXXX.5. EXTRACTION OF STATE VECTORS AT
% ANY OTHER JULIAN DATE OR TIME NOT POSSIBLE*
% The range of state vector matrix is from 1 Jan 2000 to 31 Dec 2099
% CURRENTLY SUPPORTED PLANETS:
% EARTH
% JUPITER
% SATURN
% URANUS
%% INPUTS
% State_Vector_date_and_time = Date and time at which we want the [yyyy mm dd 0 0 0]
% state vector of the desired planet. * READ ABOVE FOR EXTRACTION CONSTRAINTS*
% Planet = Desired planet for which we want to extract state vector
%% Creator:- ANKUR DEVRA 
% Develope Date - 4 October 2022
% Iteration 1 -
%% Calculations
format long
if Planet == "Earth"
    State_vectors = readmatrix("Earth_state_vectors_1_jan_2000_to_31_dec_2099.txt");
    Julian_dates = rmmissing(State_vectors(:,1)); % julian dates from 1 jan 2000 to 31 dec 2099
    X = rmmissing(State_vectors(:,3)); % km, x position 
    Y = rmmissing(State_vectors(:,4)); % km, y position 
    Z = rmmissing(State_vectors(:,5)); % km, z position 
    VX = rmmissing(State_vectors(:,6)); % km, vx velocity 
    VY = rmmissing(State_vectors(:,7)); % km, vy velocity 
    VZ = rmmissing(State_vectors(:,8)); % km, vz velocity 
    Departure_date_julian = juliandate(State_Vector_date_and_time); % departure date from earth julian date
    [row_departure,column_departure] = find(Julian_dates==Departure_date_julian); % julian date index of departure
    State_vectors_JPL_Horizons = [X(row_departure,column_departure) Y(row_departure,column_departure) Z(row_departure,column_departure) VX(row_departure,column_departure) VY(row_departure,column_departure) VZ(row_departure,column_departure)]; % km, km/s state vector
elseif Planet == "Jupiter"
    State_vectors = readmatrix("Jupiter_state_vectors_1_jan_2000_to_31_dec_2099.txt");
    Julian_dates = rmmissing(State_vectors(:,1)); % julian dates from 1 jan 2000 to 31 dec 2099
    X = rmmissing(State_vectors(:,3)); % km, x position 
    Y = rmmissing(State_vectors(:,4)); % km, y position 
    Z = rmmissing(State_vectors(:,5)); % km, z position 
    VX = rmmissing(State_vectors(:,6)); % km, vx velocity 
    VY = rmmissing(State_vectors(:,7)); % km, vy velocity 
    VZ = rmmissing(State_vectors(:,8)); % km, vz velocity 
    Departure_date_julian = juliandate(State_Vector_date_and_time); % departure date from earth julian date
    [row_departure,column_departure] = find(Julian_dates==Departure_date_julian); % julian date index of departure
    State_vectors_JPL_Horizons = [X(row_departure,column_departure) Y(row_departure,column_departure) Z(row_departure,column_departure) VX(row_departure,column_departure) VY(row_departure,column_departure) VZ(row_departure,column_departure)]; % km, km/s state vector
elseif Planet == "Saturn"
    State_vectors = readmatrix("Saturn_state_vectors_1_jan_2000_to_31_dec_2099.txt");
    Julian_dates = rmmissing(State_vectors(:,1)); % julian dates from 1 jan 2000 to 31 dec 2099
    X = rmmissing(State_vectors(:,3)); % km, x position 
    Y = rmmissing(State_vectors(:,4)); % km, y position 
    Z = rmmissing(State_vectors(:,5)); % km, z position 
    VX = rmmissing(State_vectors(:,6)); % km, vx velocity 
    VY = rmmissing(State_vectors(:,7)); % km, vy velocity 
    VZ = rmmissing(State_vectors(:,8)); % km, vz velocity 
    Departure_date_julian = juliandate(State_Vector_date_and_time); % departure date from earth julian date
    [row_departure,column_departure] = find(Julian_dates==Departure_date_julian); % julian date index of departure
    State_vectors_JPL_Horizons = [X(row_departure,column_departure) Y(row_departure,column_departure) Z(row_departure,column_departure) VX(row_departure,column_departure) VY(row_departure,column_departure) VZ(row_departure,column_departure)]; % km, km/s state vector
elseif Planet == "Uranus"
    State_vectors = readmatrix("Uranus_state_vectors_1_jan_2000_to_31_dec_2099.txt");
    Julian_dates = rmmissing(State_vectors(:,1)); % julian dates from 1 jan 2000 to 31 dec 2099
    X = rmmissing(State_vectors(:,3)); % km, x position 
    Y = rmmissing(State_vectors(:,4)); % km, y position 
    Z = rmmissing(State_vectors(:,5)); % km, z position 
    VX = rmmissing(State_vectors(:,6)); % km, vx velocity 
    VY = rmmissing(State_vectors(:,7)); % km, vy velocity 
    VZ = rmmissing(State_vectors(:,8)); % km, vz velocity 
    Departure_date_julian = juliandate(State_Vector_date_and_time); % departure date from earth julian date
    [row_departure,column_departure] = find(Julian_dates==Departure_date_julian); % julian date index of departure
    State_vectors_JPL_Horizons = [X(row_departure,column_departure) Y(row_departure,column_departure) Z(row_departure,column_departure) VX(row_departure,column_departure) VY(row_departure,column_departure) VZ(row_departure,column_departure)]; % km, km/s state vector
end
end
