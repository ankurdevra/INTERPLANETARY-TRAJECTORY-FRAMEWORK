function [Greenwich_Siderial_Time,Local_Siderial_Time,Julian_Day_Number] = UT_Time_and_Date_to_Greenwich_Siderial_Time(day_month_year,UT_hours_miniute_seconds,East_Longitude_of_a_place)
% The following code converts (day month year) and UT (hours miniute seconds) 
% To Greenwich Siderial Time and Julian day number.
% REQUIRED INPUTS:
% day_month_year = [day month year] with following acceptable range of
% values:-
% 1901=< year =<2099
%    1=< month =<12
%    1=< day =<31
% UT_hours_miniute_seconds = Universal Time in [hours minite seconds]
% East_Longitude_of_a_place = east longitude of a place on earth deg
% East longitude is measured positive eastward from greenwich meridian
% OUTPUT:
% Greenwich_Siderial_Time =  deg Greenwich Siderial time at the specified Date and UT
% Julian_Day_Number = Julian day number at desired Date and UT
% Local_Siderial_Time = local siderial time of a place on earth in deg
% (need to convert degrees to hrs); deg value is the right ascension of a
% celestial body lying on that places meridian
%% Creator:- ANKUR DEVRA 
% Develope Date - 24 March 2022
% Iteration 1 - 3 July 2022 ( added local siderial time )
%% Starting data
day = day_month_year(1); % day
month = day_month_year(2); % month
year = day_month_year(3); % year 
UT_in_hours = (UT_hours_miniute_seconds(1)) + ((UT_hours_miniute_seconds(2))/60) + ((UT_hours_miniute_seconds(3))/3600); % Universal time in hrs
Big_lambda = East_Longitude_of_a_place; % deg, east longitude of a place on earth
%% Calculations
J0 = (367*year) - fix((7/4)*(year+fix((month+9)/12))) + fix((275*month)/9) + day + 1721013.5; % Julian day number at 0 h UT
JD = J0 + UT_in_hours/24;% Julian day number at desired Date and UT
T0 = (J0-2451545)/(36525); % Time T0 in Julian centuries between Julian day J0 and J2000 epoch
theta_G0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0*T0 - (2.583*10^(-8))*T0*T0*T0;% deg Greenwich siderial time at 0 h UT
% Fix theta_G0 so that it lies in 0 to 360 range
if theta_G0 >= 360
    theta_G0 = theta_G0 - (fix(theta_G0/360))*360;
elseif theta_G0 < 0
    theta_G0 = theta_G0 - ((fix(theta_G0/360))-1)*360;
end
theta_G = theta_G0 + (360.98564726)*(UT_in_hours/24); % deg Greenwich siderial time at current UT
theta = theta_G + Big_lambda; % deg, local siderial of a place on earth
% Fix theta so that it lies in 0 to 360 range
if theta >= 360
    theta = theta - (fix(theta/360))*360;
elseif theta < 0
    theta = theta - ((fix(theta/360))-1)*360;
end
%% OUTPUT
Greenwich_Siderial_Time = theta_G;% deg Greenwich siderial time at current UT
Local_Siderial_Time = theta; % deg local siderial time of a place on earth (need to convert from deg to time)
Julian_Day_Number = JD;% Julian day number at desired Date and UT
end