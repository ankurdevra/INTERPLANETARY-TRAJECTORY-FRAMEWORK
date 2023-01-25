%% Earth Jupiter TRAJECTORY Design
clc;clear
Launch_date_and_time = [2036 8 24 0 0 0]; % [yyyy mm dd UT_hrs UT_min UT_sec]
[Earth_Departure_State_vectors_JPL_Horizons] = State_Vector_of_Planets_at_specified_Date_and_Time(Launch_date_and_time,'Earth');
Jupiter_arrival_year = 2038;
for i = 1:12
    [Arrival_date_and_time_vector(i,:)] = [Jupiter_arrival_year i 1 0 0 0];
    [Jupiter_State_vectors_JPL_Horizons(i,:)] = State_Vector_of_Planets_at_specified_Date_and_Time(Arrival_date_and_time_vector(i,:),'Jupiter');
    [elapsed_time_days(i,:)] = juliandate(Arrival_date_and_time_vector(i,:)) - juliandate(Launch_date_and_time) ;
    [Angular_momentum(i,:),Inclination(i,:),Eccentricity(i,:),RAAN(i,:),Argument_of_Perigee(i,:),True_anomaly(i,:),Radial_velocity(i,:),Time_period_of_orbit_hrs(i,:),Radius_perigee(i,:),Radius_apogee(i,:),Semimajor_axis(i,:),Semilatus_rectum(i,:),Eccentric_anomaly(i,:),Mean_anomaly(i,:),Orbits_in_a_day(i,:),state_vector_at_point_1(i,:),state_vector_at_point_2(i,:)] = Heliocentric_Orbital_elements_from_lamberts_problem(Earth_Departure_State_vectors_JPL_Horizons(1:3),Jupiter_State_vectors_JPL_Horizons(i,1:3),24*elapsed_time_days(i,:),'prograde');


end
% EARTH STATE VECTOR DATA
Earth_state_vectors = readmatrix("Earth_state_vectors_1_jan_2000_to_31_dec_2099.txt");
Jupiter_state_vectors = readmatrix("Jupiter_state_vectors_1_jan_2000_to_31_dec_2099.txt");
Julian_dates = rmmissing(Earth_state_vectors(:,1)); % julian dates from 1 sept 2036 to 31 dec 2039
X_Earth = rmmissing(Earth_state_vectors((find(Julian_dates==juliandate(Launch_date_and_time))):(365+find(Julian_dates==juliandate(Launch_date_and_time))),3)); % km, x position of Earth vector
Y_Earth = rmmissing(Earth_state_vectors((find(Julian_dates==juliandate(Launch_date_and_time))):(365+find(Julian_dates==juliandate(Launch_date_and_time))),4)); % km, y position of Earth vector
Z_Earth = rmmissing(Earth_state_vectors((find(Julian_dates==juliandate(Launch_date_and_time))):(365+find(Julian_dates==juliandate(Launch_date_and_time))),5)); % km, z position of Earth vector

% JUPITER STATE VECTOR DATA
X_Jupiter = rmmissing(Jupiter_state_vectors((find(Julian_dates==juliandate(Launch_date_and_time))):(12*365+find(Julian_dates==juliandate(Launch_date_and_time))),3)); % km, x position of Jupiter vector
Y_Jupiter = rmmissing(Jupiter_state_vectors((find(Julian_dates==juliandate(Launch_date_and_time))):(12*365+find(Julian_dates==juliandate(Launch_date_and_time))),4)); % km, y position of Jupiter vector
Z_Jupiter = rmmissing(Jupiter_state_vectors((find(Julian_dates==juliandate(Launch_date_and_time))):(12*365+find(Julian_dates==juliandate(Launch_date_and_time))),5)); % km, z position of Jupiter vector

[X,Y,Z,VX,VY,VZ] = RK4_Gravitational_numerical_integrator(state_vector_at_point_1(1,1:end),2*Time_period_of_orbit_hrs(i),'Sun');
%%
plot(X,Y);grid on;hold on;
%%
plot(X_Earth,Y_Earth);
%%
plot(X_Jupiter,Y_Jupiter);
%%
% plot3(X,Y,Z);grid on;hold on
% plot3(X_Earth,Y_Earth,Z_Earth)




%X = linspace(2038.0,2038.12,12);
% for i = 1:12
%     [Jupiter_State_vectors_JPL_Horizons(i,:)] = State_Vector_of_Planets_at_specified_Date_and_Time(Arrival_date_and_time_vector(i,:),'Jupiter');
% end
% Jupiter_State_vectors_JPL_Horizons;
% for i = 1:12
%     [elapsed_time(i,:)] = juliandate(Arrival_date_and_time_vector(i,:)) - juliandate(Launch_date_and_time) ;
% end
% elapsed_time;
% for i = 1:12
%     [Angular_momentum(i,:),Inclination(i,:),Eccentricity(i,:),RAAN(i,:),Argument_of_Perigee(i,:),True_anomaly(i,:),Radial_velocity(i,:),Time_period_of_orbit_hrs(i,:),Radius_perigee(i,:),Radius_apogee(i,:),Semimajor_axis(i,:),Semilatus_rectum(i,:),Eccentric_anomaly(i,:),Mean_anomaly(i,:),Orbits_in_a_day(i,:),state_vector_at_point_1(i,:),state_vector_at_point_2(i,:)] = Heliocentric_Orbital_elements_from_lamberts_problem(Earth_Departure_State_vectors_JPL_Horizons(1:3),Jupiter_State_vectors_JPL_Horizons(i,1:3),24*elapsed_time(i,:),'prograde')
% end
%
% figure(1)
% plot(X,Angular_momentum);grid on
% figure(2)
% plot(X,Inclination);grid on
% figure(3)
% plot(X,Eccentricity); grid on
% figure(4)
% plot(X,RAAN);grid on
% figure(5)
% plot(X,Time_period_of_orbit_hrs/2);grid on;
