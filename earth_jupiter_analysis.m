clc;clear
Years = 2038:1:2039;
Months = 1:12;
diff = (Years(end)-Years(1))+1;
X_axis_vals = linspace(Years(1),Years(end)+1,diff*length(Months));
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity(i,j),arrival_hyperbolic_excess_velocity(i,j)]  = Earth_Jupiter_Hyperbolic_Velocity_Analysis([2036 8 24 0 0 0] ,[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date = departure_hyperbolic_excess_velocity;D_Date_ = Departure_Date';D_Date_ = (D_Date_(:))';
Arrival_Date = arrival_hyperbolic_excess_velocity;A_Date_ = Arrival_Date';A_Date_ = (A_Date_(:))';
%
figure(1)
plot(X_axis_vals,D_Date_,'k');grid on;hold on
plot(X_axis_vals,A_Date_,'r')
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Sept_1(i,j),arrival_hyperbolic_excess_velocity_Date_Sept_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[2036 9 1 0 0 0],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Sept_1 = departure_hyperbolic_excess_velocity_Date_Sept_1;D_Date_Sept_1 = Departure_Date_Sept_1';D_Date_Sept_1 = (D_Date_Sept_1(:))';
Arrival_Date_Sept_1 = arrival_hyperbolic_excess_velocity_Date_Sept_1;A_Date_Sept_1 = Arrival_Date_Sept_1';A_Date_Sept_1 = (A_Date_Sept_1(:))';
figure(2)
plot(X_axis_vals,D_Date_Sept_1);hold on;grid on;
plot(X_axis_vals,A_Date_Sept_1); % case 17
%%
clc
[departure_hyperbolic_excess_velocity,arrival_hyperbolic_excess_velocity,departure_hyperbolic_velocity_vec,arrival_hyperbolic_velocity_vec] = Earth_Jupiter_Hyperbolic_Velocity_Analysis([2036 9 1 0 0 0],[2039 3 1 0 0 0]);
TOF = (juliandate([2039 3 1 0 0 0])-juliandate([2036 9 1 0 0 0]))*24; % hrs time of flight
[Earth_state_vectors_JPL_Horizons] = State_Vector_of_Planets_at_specified_Date_and_Time([2036 9 1 0 0 0],'Earth');
vec = [Earth_state_vectors_JPL_Horizons(1:3),departure_hyperbolic_velocity_vec];
[Angular_momentum,Inclination,Eccentricity,RAAN,Argument_of_Perigee,True_anomaly,Radial_velocity,Time_period_of_orbit_hrs,Radius_perigee,Radius_apogee,Semimajor_axis,Semilatus_rectum,Eccentric_anomaly,Mean_anomaly,Orbits_in_a_day] = Heliocentric_Orbital_elements_from_State_vectors([Earth_state_vectors_JPL_Horizons(1:3),departure_hyperbolic_velocity_vec]);



