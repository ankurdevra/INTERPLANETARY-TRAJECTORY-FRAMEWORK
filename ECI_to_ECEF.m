function [ECEF_State_Vector] = ECI_to_ECEF(ECI_State_Vector_Position,day_month_year,UT_hours_miniute_seconds)
% The following code converts ECI position coordinates to ECEF position coordinate at the
% specified date and UT time.
% REQUIRED INPUT:
% ECI_State_Vector_Position = [X Y Z] km current position coordinates in
% ECI frame.
% day_month_year = [day month year] with following acceptable range of
% values:-
% 1901=< year =<2099
%    1=< month =<12
%    1=< day =<31
% UT_hours_miniute_seconds = Universal Time in [hours miniute seconds]
% OUTPUT:
% ECEF_State_Vector = Converted ECEF coordinates at the specifie date and UT.
%% Calculation
% calls function UT_Time_and_Date_to_Greenwich_Siderial_Time to convert
% given date and UT to Greenwich Siderial Time.
% function [Greenwich_Siderial_Time,Julian_Day_Number] = UT_Time_and_Date_to_Greenwich_Siderial_Time(day_month_year,UT_hours_miniute_seconds)
[theta_Gt,~] = UT_Time_and_Date_to_Greenwich_Siderial_Time(day_month_year,UT_hours_miniute_seconds);
%% OUTPUT
% DCM TO CONVERT FROM ECI TO ECEF IS:-
% [EN] = | cosd(theta_Gt) sind(theta_Gt)  0 |
%        | -sind(theta_Gt) cosd(theta_Gt) 0 |
%        |       0               0        1 |
% theta_Gt = deg Greenwich Siderial time at the specified Date and UT
ECEF_State_Vector = [cosd(theta_Gt) sind(theta_Gt) 0;-sind(theta_Gt) cosd(theta_Gt) 0;0 0 1]*ECI_State_Vector_Position';
end