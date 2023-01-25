%% ANKUR DEVRA MAE 460 EARTH-JUPITER-SATURN-URANUS TRAJECTORY ANALYSIS
% The following code does trajectory analysis to launch a space probe from
% earth to the outer solar system. Gravity assist around Jupiter, Saturn and Uranus
% was analyzed and subsequent trajectory was modeled. 
% The global frame of reference for trajectory plotting and analysis was
% done with regards to ICRF frame. Furthermore, planetocentric frames of
% reference around Earth, Jupiter, Saturn and Uranus are also considered
% during the analysis. The trajectory was numerically plotted using a RK4
% gravitational numerical integrator.
clc;clear;format long
%%
% save('Final_Trajectory_Analysis')
% load('Final_Trajectory_Analysis')
%% Earth to Jupiter Analysis
Earth_Launch_Date = juliandate([2036 1 1]); % Trial Launch date of spacecraft [yyyy mm dd], 1 jan 2036
Earliest_Jupiter_Arrival_Date = juliandate([2038 1 1]); % Earliest spacecraft arrival at Jupiter [yyyy mm dd], 1 jan 2038
Latest_Jupiter_Arrival_Date = juliandate([2039 12 31]); % Latest spacecraft arrival at Jupiter [yyyy mm dd], 31 dec 2039
Launch_Date_range = juliandate([2036 1 1]):10:juliandate([2036 12 31]); % Julian day range of launch dates every 10 days apart
Arrival_Date_range = Earliest_Jupiter_Arrival_Date:10:Latest_Jupiter_Arrival_Date; % Julian day range of Jupiter arrival dates every 10 day apart
matrix_string_values = [];
matrix = [];
% This loop calculates the departure hyperbolic excess velocity and arrival
% hyperbolic excess velocity between earth and jupiter for the above
% specified launch and arrival dates.
for i = 1:length(Launch_Date_range)
    for j = 1:length(Arrival_Date_range)    
        [departure_hyperbolic_excess_velocity(i,j),arrival_hyperbolic_excess_velocity(i,j),TOF_hrs(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[datevec(datetime(Launch_Date_range(i),'convertfrom','juliandate','Format','yyyy MM dd HH mm ss'))],[datevec(datetime(Arrival_Date_range(j),'convertfrom','juliandate','Format','yyyy MM dd HH mm ss'))]);
        Launch_Dates(i) = (datetime(Launch_Date_range(i),'convertfrom','juliandate','Format','default'));
        Arrival_Dates(j) = (datetime(Arrival_Date_range(j),'convertfrom','juliandate','Format','default'));
        matrix_string_values = [matrix_string_values;string(Launch_Dates(i)),string(Arrival_Dates(j)),(departure_hyperbolic_excess_velocity(i,j)),(arrival_hyperbolic_excess_velocity(i,j)),double(TOF_hrs(i,j))];
        matrix = [matrix;Launch_Date_range(i),Arrival_Date_range(j),departure_hyperbolic_excess_velocity(i,j),arrival_hyperbolic_excess_velocity(i,j),TOF_hrs(i,j)];
    end
end
[row,~] = size(matrix_string_values);
iterations = 1:row;
Date_Velocity_table = [iterations' (matrix_string_values)];
a2t = (array2table(Date_Velocity_table,"VariableNames",["Index","Earth Launch Date","Jupiter Arrival Date","Earth Hyperbolic Departure Velocity (km/sec)","Jupiter Hyperbolic Arrival Velocity (km/sec)","Time of Flight (hrs)"]));disp(a2t); % condenses the output data into a table
Earth_Departure_Julian_Dates = matrix(:,1); % Earth departure julian date matrix
Jupiter_Arrival_Julian_Dates = matrix(:,2); % Jupiter arrival julian date matrix
Earth_Departure_Hyperbolic_excess_velocity = matrix(:,3); % km/sec, Earth departure hyperbolic excess velocity
Jupiter_Arrival_Hyperbolic_excess_velocity = matrix(:,4); % km/sec, Jupiter arrival hyperbolic excess velocity
Earth_Jupiter_ToF = matrix(:,5); % hrs, Time of flight between Earth and Jupiter
Earth_Jupiter_Velocity_matrix = [Earth_Departure_Julian_Dates,Jupiter_Arrival_Julian_Dates,Earth_Departure_Hyperbolic_excess_velocity,Jupiter_Arrival_Hyperbolic_excess_velocity,Earth_Jupiter_ToF]; % matrix of earth departure date and velocity and jupiter arrival date and velocity with time of flight
% Plots the variation of earth departure and jupiter arrival julian dates
% with respect to earth departure and jupiter arrival hyperbolic excess
% velocity.
figure(1)
plot3(Earth_Departure_Julian_Dates,Earth_Departure_Hyperbolic_excess_velocity,Jupiter_Arrival_Hyperbolic_excess_velocity,'b.'); grid on;hold on
plot3(Jupiter_Arrival_Julian_Dates,Earth_Departure_Hyperbolic_excess_velocity,Jupiter_Arrival_Hyperbolic_excess_velocity,'r.');
xlabel('Earth Departure Julian Dates/Jupiter arrival Julian Dates [x axis]');zlabel('Jupiter Arrival Hyperbolic excess velocity [km/sec] [z axis]');ylabel('Earth departure hyperbolic excess velocity [km/sec] [y axis]')
title('Earth departure and Jupiter arrival velocity analysis with 10 day intervals')
legend('Earth hyperbolic departure velocity','Jupiter hyperbolic arrival velocity ')
% This segment filter the out data with regards to chosen contraints. In
% our analysis a constraint on earth hyperbolic departure velocity was
% placed, namely that the earth hyperbolic departure velocity should be
% between 8 and 20 km/sec.
Desired_Earth_Jupiter_Date_and_Velocity_matrix = Earth_Jupiter_Velocity_matrix(Earth_Jupiter_Velocity_matrix(:,3)>=8 & Earth_Jupiter_Velocity_matrix(:,3)<=20,:); % fiters out the sub matrix containing only earth departure hyperbolic velocity less than 20 km/sec
Desired_Earth_Launch_Dates = string(datetime(Desired_Earth_Jupiter_Date_and_Velocity_matrix(:,1),'convertfrom','juliandate','Format','default'));
Desired_Jupiter_Arrival_Dates = string(datetime(Desired_Earth_Jupiter_Date_and_Velocity_matrix(:,2),'convertfrom','juliandate','Format','default'));
Desired_Earth_Jupiter_Date_and_Velocity_string_matrix = [Desired_Earth_Launch_Dates,Desired_Jupiter_Arrival_Dates,Desired_Earth_Jupiter_Date_and_Velocity_matrix(:,3),Desired_Earth_Jupiter_Date_and_Velocity_matrix(:,4),Desired_Earth_Jupiter_Date_and_Velocity_matrix(:,5)];
a2t1 = (array2table(Desired_Earth_Jupiter_Date_and_Velocity_string_matrix,"VariableNames",["Earth Launch Date","Jupiter Arrival Date","Earth Hyperbolic Departure Velocity (km/sec)","Jupiter Hyperbolic Arrival Velocity (km/sec)","Time of Flight (hrs)"])); 
a2t1 = table(a2t1,'VariableNames',{'Earth Depature and Jupiter Arrival Hyperbolic Excess Velocity'});disp(a2t1);
% Refined analyisis of optimum Earth to Jupiter launch dates
Earliest_Feasable_Earth_Launch_Date = Desired_Earth_Jupiter_Date_and_Velocity_matrix(1,1); % julian date of earliest feasable earth departure
Latest_Feasable_Earth_Launch_Date = Desired_Earth_Jupiter_Date_and_Velocity_matrix(end,1); % julian date of latest feasable earth departure
Feasable_Earth_Departure_Vector = Earliest_Feasable_Earth_Launch_Date:1:Latest_Feasable_Earth_Launch_Date; % julian date vector, one day apart, of feasable earth departure
Earliest_Feasable_Jupiter_Arrival_Date = Desired_Earth_Jupiter_Date_and_Velocity_matrix(1,2); % julian date of earliest feasable jupiter arrival
Latest_Feasable_Jupiter_Arrival_Date = Desired_Earth_Jupiter_Date_and_Velocity_matrix(end,2); % julian date of latest feasable jupiter arrival
Feasable_Jupiter_Arrival_Vector = Earliest_Jupiter_Arrival_Date:1:Latest_Jupiter_Arrival_Date;%Earliest_Feasable_Jupiter_Arrival_Date:1:Latest_Feasable_Jupiter_Arrival_Date; % julian date vector, one day apart, of feasable jupiter arrival
% This loop calculates the departure hyperbolic excess velocity and arrival
% hyperbolic excess velocity between earth and jupiter for the above
% specified refined launch and arrival date ranges.
Fesable_matrix_string_values = [];
Fesable_matrix = [];
for i = 1:length(Feasable_Earth_Departure_Vector)
    for j = 1:length(Feasable_Jupiter_Arrival_Vector)
        [Fesable_departure_hyperbolic_excess_velocity(i,j),Fesable_arrival_hyperbolic_excess_velocity(i,j),Fesable_TOF_hrs(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[datevec(datetime(Feasable_Earth_Departure_Vector(i),'convertfrom','juliandate','Format','yyyy MM dd HH mm ss'))],[datevec(datetime(Feasable_Jupiter_Arrival_Vector(j),'convertfrom','juliandate','Format','yyyy MM dd HH mm ss'))]);
        Fesable_Launch_Dates(i) = (datetime(Feasable_Earth_Departure_Vector(i),'convertfrom','juliandate','Format','default'));
        Fesable_Arrival_Dates(j) = (datetime(Feasable_Jupiter_Arrival_Vector(j),'convertfrom','juliandate','Format','default'));
        Fesable_matrix_string_values = [Fesable_matrix_string_values;string(Fesable_Launch_Dates(i)),string(Fesable_Arrival_Dates(j)),(Fesable_departure_hyperbolic_excess_velocity(i,j)),(Fesable_arrival_hyperbolic_excess_velocity(i,j)),double(Fesable_TOF_hrs(i,j))];
        Fesable_matrix = [Fesable_matrix;Feasable_Earth_Departure_Vector(i),Feasable_Jupiter_Arrival_Vector(j),Fesable_departure_hyperbolic_excess_velocity(i,j),Fesable_arrival_hyperbolic_excess_velocity(i,j),Fesable_TOF_hrs(i,j)];
    end
end
[Fesable_row,~] = size(Fesable_matrix_string_values);
Fesable_iterations = 1:Fesable_row;
Fesable_Date_Velocity_table = [Fesable_iterations' (Fesable_matrix_string_values)];
a2t2 = (array2table(Fesable_Date_Velocity_table,"VariableNames",["Index","Earth Launch Date","Jupiter Arrival Date","Earth Hyperbolic Departure Velocity (km/sec)","Jupiter Hyperbolic Arrival Velocity (km/sec)","Time of Flight (hrs)"]));
a2t2 = table(a2t2,'VariableNames',{'Refined Earth Depature and Jupiter Arrival Hyperbolic Velocity'});disp(a2t2);
% Plots the variation of earth departure and jupiter arrival julian dates
% with respect to earth departure and jupiter arrival hyperbolic excess
% velocity using refined date ranges.
figure(2)
plot3(Feasable_Earth_Departure_Vector,Fesable_departure_hyperbolic_excess_velocity,Fesable_arrival_hyperbolic_excess_velocity,'b.'); grid on;hold on
plot3(Feasable_Jupiter_Arrival_Vector,Fesable_departure_hyperbolic_excess_velocity,Fesable_arrival_hyperbolic_excess_velocity,'r.');
xlabel('Earth Departure Julian Dates/Jupiter arrival Julian Dates [x axis]');zlabel('Jupiter Arrival Hyperbolic excess velocity [km/sec] [z axis]');ylabel('Earth departure hyperbolic excess velocity [km/sec] [y axis]')
title(' Refined Earth departure and Jupiter arrival velocity analysis')
legend('Earth hyperbolic departure velocity','Jupiter hyperbolic arrival velocity')
% ORBIT PLOTTING OF EARTH AND JUPITER
Orbit_plotting_vector_julian_dates = Earth_Launch_Date:1:Latest_Jupiter_Arrival_Date; % julian date vector for orbit plotting
for i = 1:length(Orbit_plotting_vector_julian_dates)
    [Earth_Position_ICRF(i,:),Earth_Velocity_ICRF(i,:)] = planetEphemeris(Orbit_plotting_vector_julian_dates(i),'Sun','Earth','432t');
    Earth_Position_ICRF_norm(i,:) = norm(Earth_Position_ICRF(i,:));
    Earth_Velocity_ICRF_norm(i,:) = norm(Earth_Velocity_ICRF(i,:));
    [Jupiter_Position_ICRF(i,:),Jupiter_Velocity_ICRF(i,:)] = planetEphemeris(Orbit_plotting_vector_julian_dates(i),'Sun','Jupiter','432t');
    Jupiter_Position_ICRF_norm(i,:) = norm(Jupiter_Position_ICRF(i,:));
    Jupiter_Velocity_ICRF_norm(i,:) = norm(Jupiter_Velocity_ICRF(i,:));
end
figure(3)
plot3(Earth_Position_ICRF(:,1),Earth_Position_ICRF(:,2),Earth_Position_ICRF(:,3),'b'); % Earth orbit in ICRF frame
grid on; hold on;
plot3(Jupiter_Position_ICRF(:,1),Jupiter_Position_ICRF(:,2),Jupiter_Position_ICRF(:,3),'r'); % Jupiter orbit in ICRF frame
plot3(0,0,0, 'o', 'MarkerFaceColor', 'y') % yellow dot representing sun, not to scale
xlabel('X [km]');ylabel('Y [km]');zlabel('Z [km]')
% Plotting earth jupiter trajectory
Achivable_Earth_Departure_Hyperbolic_Velocity = 16; % km/sec, Achivable Earth Departure Hyperbolic Velocity,to be set by us
[~,~,index] = unique((abs(Fesable_matrix(:,3)-Achivable_Earth_Departure_Hyperbolic_Velocity))); % achievable launch date based on our achivable velocity
Achivable_Earth_Departure_Hyperbolic_Velocity_Date_vector = Fesable_matrix(index==1,:); % vector of date and velocity required at achievable launch date
Achievable_Earth_Launch_Date = Achivable_Earth_Departure_Hyperbolic_Velocity_Date_vector(1);
Achievable_Jupiter_Arrival_Date = Achivable_Earth_Departure_Hyperbolic_Velocity_Date_vector(2);
Achievable_Earth_Launch_Date_string = datetime(Achivable_Earth_Departure_Hyperbolic_Velocity_Date_vector(1),'convertfrom','juliandate','Format','default');
Achievable_Jupiter_Arrival_Date_string = datetime(Achivable_Earth_Departure_Hyperbolic_Velocity_Date_vector(2),'convertfrom','juliandate','Format','default');
[Position_Earth_Launch,Velocity_Earth_Launch] = planetEphemeris(Achievable_Earth_Launch_Date,'Sun','Earth','432t'); % km/sec, earth position and velocity at time of launch 
[Position_Jupiter_Arrival,Velocity_Jupiter_Arrival] = planetEphemeris(Achievable_Jupiter_Arrival_Date,'Sun','Jupiter','432t');% km/sec, jupiter position and velocity at time of launch 
[Angular_momentum_Earth_Jupiter_spacecraft_orbit,~,Eccentricity_Earth_Jupiter_spacecraft_orbit,~,~,True_anomaly_Earth_Jupiter_spacecraft_orbit,~,~,~,~,Semimajor_axis_Earth_Jupiter_spacecraft_orbit,~,~,~,~,state_vector_at_point_1_Earth_Jupiter_spacecraft_orbit,state_vector_at_point_2_Earth_Jupiter_spacecraft_orbit] = Heliocentric_Orbital_elements_from_lamberts_problem(Position_Earth_Launch,Position_Jupiter_Arrival,(Achievable_Jupiter_Arrival_Date-Achievable_Earth_Launch_Date)*24,'prograde')
[X_Earth_Jupiter_spacecraft_orbit,Y_Earth_Jupiter_spacecraft_orbit,Z_Earth_Jupiter_spacecraft_orbit,VX_Earth_Jupiter_spacecraft_orbit,VY_Earth_Jupiter_spacecraft_orbit,VZ_Earth_Jupiter_spacecraft_orbit] = RK4_Gravitational_numerical_integrator(state_vector_at_point_1_Earth_Jupiter_spacecraft_orbit,(Achievable_Jupiter_Arrival_Date-Achievable_Earth_Launch_Date)*24,'Sun');
plot3(Position_Earth_Launch(1),Position_Earth_Launch(2),Position_Earth_Launch(3), 'o', 'MarkerFaceColor', 'b');
plot3(Position_Jupiter_Arrival(1),Position_Jupiter_Arrival(2),Position_Jupiter_Arrival(3), 'o', 'MarkerFaceColor', 'r')
plot3(X_Earth_Jupiter_spacecraft_orbit,Y_Earth_Jupiter_spacecraft_orbit,Z_Earth_Jupiter_spacecraft_orbit,'k')
legend('Earth Orbit','Jupiter orbit','Sun','Earth Position at Launch','Jupiter Position at Arrival','Earth-Jupiter Spacecraft Coasting Orbit')
% Jupiter Gravity assist analysis
mu_Sun = 132712440041.93938; % km^3/s^2 sun gravitational parameter
mu_Jupiter = 126686531.900; % km^3/sec^2 gravitational prameter of jupiter
h_1_Jupiter = Angular_momentum_Earth_Jupiter_spacecraft_orbit; % km^2/sec, Angular momentum of heliocentric approach trajectory
e_1_Jupiter = Eccentricity_Earth_Jupiter_spacecraft_orbit; % eccentricity of heliocentric approach trajectory
theta_1_Jupiter = True_anomaly_Earth_Jupiter_spacecraft_orbit; % deg, True anomaly of heliocentric approach trajectory
V_Transverse_1_Spacecraft_Jupiter = (mu_Sun/h_1_Jupiter)*(1+(e_1_Jupiter*cosd(theta_1_Jupiter))); % km/sec, transverse component of spacecraft velocity at inbound crossing
V_Radial_1_Spacecraft_Jupiter = (mu_Sun/h_1_Jupiter)*((e_1_Jupiter*sind(theta_1_Jupiter))); % km/sec, radial component of spacecraft velocity at inbound crossing
V_1_spacecraft_Jupiter_vector = [V_Transverse_1_Spacecraft_Jupiter V_Radial_1_Spacecraft_Jupiter]; % km/sec, [u_v u_s] heliocentric velocity of spacecraft at inbound crossing
V_Jupiter_vector = [norm(Velocity_Jupiter_Arrival) 0]; % km/sec, [u_v u_s] heliocentric velocity of Jupiter
v_infinity_1_Jupiter_vector = V_1_spacecraft_Jupiter_vector-V_Jupiter_vector; % km/sec,[u_v u_s] hyperbolic excess velocity of spacecraft at inbound crossing
v_infinity__1_Jupiter = norm(v_infinity_1_Jupiter_vector); % km/sec, inbound hyperbolic excess velocity of the spacecraft realtive to Jupiter
r_p_Jupiter = 300000; % km, periapsis radius of spacecraft around jupiter during flyby from center of jupiter
h_Jupiter = r_p_Jupiter.*sqrt(((v_infinity__1_Jupiter)^2)+((2*mu_Jupiter)./(r_p_Jupiter))); % km^2/sec, angular momentum of planetocentric hyperbola
e_Jupiter = 1 + ((r_p_Jupiter.*((v_infinity__1_Jupiter)^2))/(mu_Jupiter)); % eccentricity of planetocentric hyperbola
delta_Jupiter = 2*asind(1./e_Jupiter); % deg, turn angle of spacecraft in flyby wrt Jupiter
theta_infinity_Jupiter = acosd(-1./e_Jupiter); % deg, true anomaly of asymptote of spacecraft in flyby wrt Jupiter
phi_1_Jupiter = atand(v_infinity_1_Jupiter_vector(2)/v_infinity_1_Jupiter_vector(1)); % deg, angle betweenn v_infinity_1 and V (Jupiter velocity vector) at inbound crossing
phi_2_Jupiter = phi_1_Jupiter-delta_Jupiter; % deg, angle betweenn v_infinity_1 and V (Jupiter velocity vector) at outbound crossing
v_infinity__2_Jupiter_vector = v_infinity__1_Jupiter.*[cosd(phi_2_Jupiter) sind(phi_2_Jupiter)]; % km/sec, [u_v u_s] hyperbolic excess velocity of spacecraft at outbound crossing
v_infinity__2_Jupiter = norm(v_infinity__2_Jupiter_vector); % km/sec, outbound hyperbolic excess velocity of the spacecraft realtive to Jupiter
V_2_spacecraft_Jupiter_vector = v_infinity__2_Jupiter_vector+V_Jupiter_vector; % km/sec, [u_v u_s] heliocentric velocity of spacecraft at outbound crossing
V_inbound_Jupiter = norm(V_1_spacecraft_Jupiter_vector); % km/sec, inbound heliocentric velocity of spacecraft wrt Jupiter
V_outbound_Jupiter = norm(V_2_spacecraft_Jupiter_vector); % km/sec, outbound heliocentric velocity of spacecraft wrt Jupiter
delta_v_Jupiter_flyby = V_outbound_Jupiter-V_inbound_Jupiter; % km/sec, delta-v change in heliocentric velocity of spacecraft afgter jupiter flyby
% Jupiter-Saturn Trajectory Analysis
Jupiter_Departure_Julian_Date = Achievable_Jupiter_Arrival_Date; % Julian date when spacecraft leaves Jupiter outbound for Saturn
Jupiter_Departure_Date = datetime(Jupiter_Departure_Julian_Date,'convertfrom','juliandate','Format','default'); % Date when spacecraft leaves Jupiter outbound for Saturn
Saturn_Arrival_Julian_Date_vector = (Jupiter_Departure_Julian_Date+100):(Jupiter_Departure_Julian_Date+1000); % Julian dates of Saturn arrival between +100 and +1000 days after Jupiter departure
for i = 1:length(Saturn_Arrival_Julian_Date_vector)
    [Saturn_Arrival_Position_Vector_ICRF(i,:),Saturn_Arrival_Velocity_Vector_ICRF(i,:)] = planetEphemeris(Saturn_Arrival_Julian_Date_vector(i),'Sun','Saturn','432t');
end
Feasable_Matrix_Jupiter_Saturn_Trajectory = [];
Feasable_Matrix_Jupiter_Saturn_Trajectory_String = [];
% Calulates the orbital parameters of the Jupiter-Saturn trajectory for a
% range of Saturn arrival dates
for i = 1:length(Saturn_Arrival_Julian_Date_vector)
    [Angular_momentum_Jupiter_Saturn_spacecraft_orbit(i),Inclination_Jupiter_Saturn_spacecraft_orbit(i),Eccentricity_Jupiter_Saturn_spacecraft_orbit(i),~,~,True_anomaly_Jupiter_Saturn_spacecraft_orbit(i),~,~,~,~,Semimajor_axis_Jupiter_Saturn_spacecraft_orbit(i),~,~,~,~,state_vector_at_point_1_Jupiter_Saturn_spacecraft_orbit(i,:),state_vector_at_point_2_Jupiter_Saturn_spacecraft_orbit(i,:)] = Heliocentric_Orbital_elements_from_lamberts_problem(Position_Jupiter_Arrival,Saturn_Arrival_Position_Vector_ICRF(i,:),(Saturn_Arrival_Julian_Date_vector(i)-Jupiter_Departure_Julian_Date)*24,'prograde');
    Feasable_Saturn_Arrival_Dates(i) = (datetime(Saturn_Arrival_Julian_Date_vector(i),'convertfrom','juliandate','Format','default'));
    Feasable_Matrix_Jupiter_Saturn_Trajectory_String = [Feasable_Matrix_Jupiter_Saturn_Trajectory_String;Angular_momentum_Jupiter_Saturn_spacecraft_orbit(i),Inclination_Jupiter_Saturn_spacecraft_orbit(i),Eccentricity_Jupiter_Saturn_spacecraft_orbit(i),True_anomaly_Jupiter_Saturn_spacecraft_orbit(i),Semimajor_axis_Jupiter_Saturn_spacecraft_orbit(i),string(Feasable_Saturn_Arrival_Dates(i)),state_vector_at_point_1_Jupiter_Saturn_spacecraft_orbit(i,:),state_vector_at_point_2_Jupiter_Saturn_spacecraft_orbit(i,:)];
    Feasable_Matrix_Jupiter_Saturn_Trajectory = [Feasable_Matrix_Jupiter_Saturn_Trajectory;Angular_momentum_Jupiter_Saturn_spacecraft_orbit(i),Inclination_Jupiter_Saturn_spacecraft_orbit(i),Eccentricity_Jupiter_Saturn_spacecraft_orbit(i),True_anomaly_Jupiter_Saturn_spacecraft_orbit(i),Semimajor_axis_Jupiter_Saturn_spacecraft_orbit(i),Saturn_Arrival_Julian_Date_vector(i),state_vector_at_point_1_Jupiter_Saturn_spacecraft_orbit(i,:),state_vector_at_point_2_Jupiter_Saturn_spacecraft_orbit(i,:)];
end
% Jupiter Uranus Trajectory analysis
Uranus_Arrival_Julian_Date_vector = (Jupiter_Departure_Julian_Date+1500):(Jupiter_Departure_Julian_Date+2500); % Julian dates of Uranus arrival between +1500 and +2500 days after Jupiter departure
for i = 1:length(Uranus_Arrival_Julian_Date_vector)
    [Uranus_Arrival_Position_Vector_ICRF(i,:),Uranus_Arrival_Velocity_Vector_ICRF(i,:)] = planetEphemeris(Uranus_Arrival_Julian_Date_vector(i),'Sun','Uranus','432t');
end
Feasable_Matrix_Jupiter_Uranus_Trajectory = [];
Feasable_Matrix_Jupiter_Uranus_Trajectory_String = [];
% Calulates the orbital parameters of the Jupiter-Uranus trajectory for a
% range of Uranus arrival dates
for i = 1:length(Uranus_Arrival_Julian_Date_vector)
    [Angular_momentum_Jupiter_Uranus_spacecraft_orbit(i),Inclination_Jupiter_Uranus_spacecraft_orbit(i),Eccentricity_Jupiter_Uranus_spacecraft_orbit(i),~,~,True_anomaly_Jupiter_Uranus_spacecraft_orbit(i),~,~,~,~,Semimajor_axis_Jupiter_Uranus_spacecraft_orbit(i),~,~,~,~,state_vector_at_point_1_Jupiter_Uranus_spacecraft_orbit(i,:),state_vector_at_point_2_Jupiter_Uranus_spacecraft_orbit(i,:)] = Heliocentric_Orbital_elements_from_lamberts_problem(Position_Jupiter_Arrival,Uranus_Arrival_Position_Vector_ICRF(i,:),(Uranus_Arrival_Julian_Date_vector(i)-Jupiter_Departure_Julian_Date)*24,'prograde');
    Feasable_Uranus_Arrival_Dates(i) = (datetime(Uranus_Arrival_Julian_Date_vector(i),'convertfrom','juliandate','Format','default'));
    Feasable_Matrix_Jupiter_Uranus_Trajectory_String = [Feasable_Matrix_Jupiter_Uranus_Trajectory_String;Angular_momentum_Jupiter_Uranus_spacecraft_orbit(i),Inclination_Jupiter_Uranus_spacecraft_orbit(i),Eccentricity_Jupiter_Uranus_spacecraft_orbit(i),True_anomaly_Jupiter_Uranus_spacecraft_orbit(i),Semimajor_axis_Jupiter_Uranus_spacecraft_orbit(i),string(Feasable_Uranus_Arrival_Dates(i)),state_vector_at_point_1_Jupiter_Uranus_spacecraft_orbit(i,:),state_vector_at_point_2_Jupiter_Uranus_spacecraft_orbit(i,:)];
    Feasable_Matrix_Jupiter_Uranus_Trajectory = [Feasable_Matrix_Jupiter_Uranus_Trajectory;Angular_momentum_Jupiter_Uranus_spacecraft_orbit(i),Inclination_Jupiter_Uranus_spacecraft_orbit(i),Eccentricity_Jupiter_Uranus_spacecraft_orbit(i),True_anomaly_Jupiter_Uranus_spacecraft_orbit(i),Semimajor_axis_Jupiter_Uranus_spacecraft_orbit(i),Uranus_Arrival_Julian_Date_vector(i),state_vector_at_point_1_Jupiter_Uranus_spacecraft_orbit(i,:),state_vector_at_point_2_Jupiter_Uranus_spacecraft_orbit(i,:)];
end
[~,~,index_] = unique((abs(Jupiter_Departure_Velocity-V_outbound_Jupiter))); % filters and Finds the departure velocity from jupiter to saturn
% optimum Jupiter saturn trajectory
Jupiter_Saturn_Trajectory_Optimum = Feasable_Matrix_Jupiter_Saturn_Trajectory(index_==1,:);
Saturn_Arrival_Date_String = datetime(Jupiter_Saturn_Trajectory_Optimum(6),'convertfrom','juliandate','Format','default'); % date when we arrive at saturn
Jupiter_Saturn_Departure_State_Vector = (Feasable_Matrix_Jupiter_Saturn_Trajectory(index_==1,7:12));
Jupiter_Saturn_Arrival_State_Vector = Feasable_Matrix_Jupiter_Saturn_Trajectory(index_==1,13:18);
[Saturn_Position_Arrival,Saturn_Velocity_Arrival] = planetEphemeris(Jupiter_Saturn_Trajectory_Optimum(6),'Sun','Saturn','432t');
[X_Jupiter_Saturn_spacecraft_orbit,Y_Jupiter_Saturn_spacecraft_orbit,Z_Jupiter_Saturn_spacecraft_orbit,VX_Jupiter_Saturn_spacecraft_orbit,VY_Jupiter_Saturn_spacecraft_orbit,VZ_Jupiter_Saturn_spacecraft_orbit] = RK4_Gravitational_numerical_integrator(Jupiter_Saturn_Departure_State_Vector,(Jupiter_Saturn_Trajectory_Optimum(6)-Jupiter_Departure_Julian_Date)*24,'Sun');
% Optimum Jupiter uranus Trajectory
Jupiter_Uranus_Trajectory_Optimum = Feasable_Matrix_Jupiter_Uranus_Trajectory(index_==1,:);
Uranus_Arrival_Date_String = datetime(Jupiter_Uranus_Trajectory_Optimum(6),'convertfrom','juliandate','Format','default'); % date when we arrive at uranus
Jupiter_Uranus_Departure_State_Vector = (Feasable_Matrix_Jupiter_Uranus_Trajectory(index_==1,7:12));
Jupiter_Uranus_Arrival_State_Vector = Feasable_Matrix_Jupiter_Uranus_Trajectory(index_==1,13:18);
[Uranus_Position_Arrival,Uranus_Velocity_Arrival] = planetEphemeris(Jupiter_Uranus_Trajectory_Optimum(6),'Sun','Uranus','432t');
[X_Jupiter_Uranus_spacecraft_orbit,Y_Jupiter_Uranus_spacecraft_orbit,Z_Jupiter_Uranus_spacecraft_orbit,VX_Jupiter_Uranus_spacecraft_orbit,VY_Jupiter_Uranus_spacecraft_orbit,VZ_Jupiter_Uranus_spacecraft_orbit] = RK4_Gravitational_numerical_integrator(Jupiter_Uranus_Departure_State_Vector,(Jupiter_Uranus_Trajectory_Optimum(6)-Jupiter_Departure_Julian_Date)*24,'Sun');
% SATURN GRAVITY ASSIST
mu_Saturn = 37931187; % km^3/sec^2 gravitational prameter of saturn
h_1_Saturn = Jupiter_Saturn_Trajectory_Optimum(1); % km^2/sec, Angular momentum of heliocentric approach trajectory to saturn
e_1_Saturn = Jupiter_Saturn_Trajectory_Optimum(3); % eccentricity of heliocentric approach trajectory
theta_1_Saturn = Jupiter_Saturn_Trajectory_Optimum(4); % deg, True anomaly of heliocentric approach trajectory
V_Transverse_1_Spacecraft_Saturn = (mu_Sun/h_1_Saturn)*(1+(e_1_Saturn*cosd(theta_1_Saturn))); % km/sec, transverse component of spacecraft velocity at inbound crossing
V_Radial_1_Spacecraft_Saturn = (mu_Sun/h_1_Saturn)*((e_1_Saturn*sind(theta_1_Saturn))); % km/sec, radial component of spacecraft velocity at inbound crossing
V_1_spacecraft_Saturn_vector = [V_Transverse_1_Spacecraft_Saturn V_Radial_1_Spacecraft_Saturn]; % km/sec, [u_v u_s] heliocentric velocity of spacecraft at inbound crossing
V_Saturn_vector = [norm(Saturn_Velocity_Arrival) 0]; % km/sec, [u_v u_s] heliocentric velocity of Jupiter
v_infinity_1_Saturn_vector = V_1_spacecraft_Saturn_vector-V_Saturn_vector; % km/sec,[u_v u_s] hyperbolic excess velocity of spacecraft at inbound crossing
v_infinity__1_Saturn = norm(v_infinity_1_Saturn_vector); % km/sec, inbound hyperbolic excess velocity of the spacecraft realtive to Saturn
r_p_Saturn = 200000; % km, periapsis radius of spacecraft around Saturn during flyby from center of Saturn
h_Saturn = r_p_Saturn.*sqrt(((v_infinity__1_Saturn)^2)+((2*mu_Saturn)./(r_p_Saturn))); % km^2/sec, angular momentum of planetocentric hyperbola
e_Saturn = 1 + ((r_p_Saturn.*((v_infinity__1_Saturn)^2))/(mu_Saturn)); % eccentricity of planetocentric hyperbola
delta_Saturn = 2*asind(1./e_Saturn); % deg, turn angle of spacecraft in flyby wrt Saturn
theta_infinity_Saturn = acosd(-1./e_Saturn); % deg, true anomaly of asymptote of spacecraft in flyby wrt uranus
phi_1_Saturn = atand(v_infinity_1_Saturn_vector(2)/v_infinity_1_Saturn_vector(1)); % deg, angle betweenn v_infinity_1 and V (Saturn velocity vector) at inbound crossing
phi_2_Saturn = phi_1_Saturn-delta_Saturn; % deg, angle betweenn v_infinity_1 and V (Saturn velocity vector) at outbound crossing
v_infinity__2_Saturn_vector = v_infinity__1_Saturn.*[cosd(phi_2_Saturn) sind(phi_2_Saturn)]; % km/sec, [u_v u_s] hyperbolic excess velocity of spacecraft at outbound crossing
v_infinity__2_Saturn = norm(v_infinity__2_Saturn_vector); % km/sec, outbound hyperbolic excess velocity of the spacecraft realtive to Saturn
V_2_spacecraft_Saturn_vector = v_infinity__2_Saturn_vector+V_Saturn_vector; % km/sec, [u_v u_s] heliocentric velocity of spacecraft at outbound crossing
V_inbound_Saturn = norm(V_1_spacecraft_Saturn_vector); % km/sec, inbound heliocentric velocity of spacecraft wrt saturn
V_outbound_Saturn = norm(V_2_spacecraft_Saturn_vector); % km/sec, outbound heliocentric velocity of spacecraft wrt saturn
delta_v_Saturn_flyby = V_outbound_Saturn-V_inbound_Saturn; % km/sec, delta-v change in heliocentric velocity of spacecraft after saturn flyby
% Uranus Gravity assist
mu_Uranus = 5793939; % km^3/sec^2 gravitational prameter of Uranus
h_1_Uranus = Jupiter_Uranus_Trajectory_Optimum(1); % km^2/sec, Angular momentum of heliocentric approach trajectory to uranus
e_1_Uranus = Jupiter_Uranus_Trajectory_Optimum(3); % eccentricity of heliocentric approach trajectory
theta_1_Uranus = Jupiter_Uranus_Trajectory_Optimum(4); % deg, True anomaly of heliocentric approach trajectory
V_Transverse_1_Spacecraft_Uranus = (mu_Sun/h_1_Uranus)*(1+(e_1_Uranus*cosd(theta_1_Uranus))); % km/sec, transverse component of spacecraft velocity at inbound crossing
V_Radial_1_Spacecraft_Uranus = (mu_Sun/h_1_Uranus)*((e_1_Uranus*sind(theta_1_Uranus))); % km/sec, radial component of spacecraft velocity at inbound crossing
V_1_spacecraft_Uranus_vector = [V_Transverse_1_Spacecraft_Uranus V_Radial_1_Spacecraft_Uranus]; % km/sec, [u_v u_s] heliocentric velocity of spacecraft at inbound crossing
V_Uranus_vector = [norm(Uranus_Velocity_Arrival) 0]; % km/sec, [u_v u_s] heliocentric velocity of Uranus
v_infinity_1_Uranus_vector = V_1_spacecraft_Uranus_vector-V_Uranus_vector; % km/sec,[u_v u_s] hyperbolic excess velocity of spacecraft at inbound crossing
v_infinity__1_Uranus = norm(v_infinity_1_Uranus_vector); % km/sec, inbound hyperbolic excess velocity of the spacecraft realtive to Uranus
r_p_Uranus = 100000; % km, periapsis radius of spacecraft around Saturn during flyby from center of Uranus
h_Uranus = r_p_Uranus.*sqrt(((v_infinity__1_Uranus)^2)+((2*mu_Uranus)./(r_p_Uranus))); % km^2/sec, angular momentum of planetocentric hyperbola
e_Uranus = 1 + ((r_p_Uranus.*((v_infinity__1_Uranus)^2))/(mu_Uranus)); % eccentricity of planetocentric hyperbola
delta_Uranus = 2*asind(1./e_Uranus); % deg, turn angle of spacecraft in flyby wrt Uranus
theta_infinity_Uranus = acosd(-1./e_Uranus); % deg, true anomaly of asymptote of spacecraft in flyby wrt Uranus
phi_1_Uranus = atand(v_infinity_1_Uranus_vector(2)/v_infinity_1_Uranus_vector(1)); % deg, angle betweenn v_infinity_1 and V (Uranus velocity vector) at inbound crossing
phi_2_Uranus = phi_1_Uranus-delta_Uranus; % deg, angle betweenn v_infinity_1 and V (Uranus velocity vector) at outbound crossing
v_infinity__2_Uranus_vector = v_infinity__1_Uranus.*[cosd(phi_2_Uranus) sind(phi_2_Uranus)]; % km/sec, [u_v u_s] hyperbolic excess velocity of spacecraft at outbound crossing
v_infinity__2_Uranus = norm(v_infinity__2_Uranus_vector); % km/sec, outbound hyperbolic excess velocity of the spacecraft realtive to Uranus
V_2_spacecraft_Uranus_vector = v_infinity__2_Uranus_vector+V_Uranus_vector; % km/sec, [u_v u_s] heliocentric velocity of spacecraft at outbound crossing
V_inbound_Uranus = norm(V_1_spacecraft_Uranus_vector); % km/sec, inbound heliocentric velocity of spacecraft wrt Uranus
V_outbound_Uranus = norm(V_2_spacecraft_Uranus_vector); % km/sec, outbound heliocentric velocity of spacecraft wrt Uranus
delta_v_Uranus_flyby = V_outbound_Uranus-V_inbound_Uranus; % km/sec, delta-v change in heliocentric velocity of spacecraft after Uranus flyby
% Plots the combined trajectory of probe 1 and 2.
figure(4)
plot3(Earth_Position_ICRF(:,1),Earth_Position_ICRF(:,2),Earth_Position_ICRF(:,3),'b');grid on; hold on; % Earth orbit in ICRF frame
plot3(Jupiter_Position_ICRF(:,1),Jupiter_Position_ICRF(:,2),Jupiter_Position_ICRF(:,3),'r'); % Jupiter orbit in ICRF frame
plot3(Saturn_Arrival_Position_Vector_ICRF(:,1),Saturn_Arrival_Position_Vector_ICRF(:,2),Saturn_Arrival_Position_Vector_ICRF(:,3),'m'); % Saturn orbit in ICRF frame
plot3(Uranus_Arrival_Position_Vector_ICRF(:,1),Uranus_Arrival_Position_Vector_ICRF(:,2),Uranus_Arrival_Position_Vector_ICRF(:,3),'b'); % Uranus orbit in ICRF frame
plot3(X_Earth_Jupiter_spacecraft_orbit,Y_Earth_Jupiter_spacecraft_orbit,Z_Earth_Jupiter_spacecraft_orbit,'k'); % Earth Jupiter Spacecraft orbit
plot3(X_Jupiter_Saturn_spacecraft_orbit,Y_Jupiter_Saturn_spacecraft_orbit,Z_Jupiter_Saturn_spacecraft_orbit,'g') % Jupiter Saturn Spcecraft orbit
plot3(X_Jupiter_Uranus_spacecraft_orbit,Y_Jupiter_Uranus_spacecraft_orbit,Z_Jupiter_Uranus_spacecraft_orbit,'b') % Jupiter Uranus Spacecraft orbit
plot3(0,0,0, 'o', 'MarkerFaceColor', 'y') % yellow dot representing sun, not to scale
plot3(Position_Earth_Launch(1),Position_Earth_Launch(2),Position_Earth_Launch(3), 'o', 'MarkerFaceColor', 'b'); % earth position at launch
plot3(Position_Jupiter_Arrival(1),Position_Jupiter_Arrival(2),Position_Jupiter_Arrival(3), 'o', 'MarkerFaceColor', 'r') % jupiter position at arrival
plot3(Saturn_Position_Arrival(1),Saturn_Position_Arrival(2),Saturn_Position_Arrival(3),'o', 'MarkerFaceColor', 'm') % saturn position at arrival
plot3(Uranus_Position_Arrival(1),Uranus_Position_Arrival(2),Uranus_Position_Arrival(3),'o','MarkerFaceColor', 'b')% Uranus position at arrival
legend('Earth Orbit','Jupiter orbit','Saturn Orbit','Uranus Orbit','Earth-Jupiter Spacecraft Coasting Orbit','Jupiter-Saturn Spacecraft Coasting Orbit','Jupiter-Uranus Spacecraft Coasting Orbit','Sun','Earth Position at Launch','Jupiter Position at Arrival','Saturn Position at Arrival','Uranus Position at Arrival')
title('Probe 1:Earth-Jupiter-Saturn Trajectory')
xlabel('X [km]');ylabel('Y [km]');zlabel('Z [km]')
% Plots the trajectory of probe 1
figure(5)
plot3(Earth_Position_ICRF(:,1),Earth_Position_ICRF(:,2),Earth_Position_ICRF(:,3),'b');grid on; hold on; % Earth orbit in ICRF frame
plot3(Jupiter_Position_ICRF(:,1),Jupiter_Position_ICRF(:,2),Jupiter_Position_ICRF(:,3),'r'); % Jupiter orbit in ICRF frame
plot3(Saturn_Arrival_Position_Vector_ICRF(:,1),Saturn_Arrival_Position_Vector_ICRF(:,2),Saturn_Arrival_Position_Vector_ICRF(:,3),'m'); % Saturn orbit in ICRF frame
plot3(X_Earth_Jupiter_spacecraft_orbit,Y_Earth_Jupiter_spacecraft_orbit,Z_Earth_Jupiter_spacecraft_orbit,'k'); % Earth Jupiter Spacecraft orbit
plot3(X_Jupiter_Saturn_spacecraft_orbit,Y_Jupiter_Saturn_spacecraft_orbit,Z_Jupiter_Saturn_spacecraft_orbit,'g') % Jupiter Saturn Spcecraft orbit
plot3(0,0,0, 'o', 'MarkerFaceColor', 'y') % yellow dot representing sun, not to scale
plot3(Position_Earth_Launch(1),Position_Earth_Launch(2),Position_Earth_Launch(3), 'o', 'MarkerFaceColor', 'b'); % earth position at launch
plot3(Position_Jupiter_Arrival(1),Position_Jupiter_Arrival(2),Position_Jupiter_Arrival(3), 'o', 'MarkerFaceColor', 'r') % jupiter position at arrival
plot3(Saturn_Position_Arrival(1),Saturn_Position_Arrival(2),Saturn_Position_Arrival(3),'o', 'MarkerFaceColor', 'm') % saturn position at arrival
title('Probe 1:Earth-Jupiter-Saturn Trajectory in ICRF frame')
xlabel('X [km]');ylabel('Y [km]');zlabel('Z [km]')
legend('Earth Orbit','Jupiter orbit','Saturn Orbit','Earth-Jupiter Spacecraft Coasting Orbit','Jupiter-Saturn Spacecraft Coasting Orbit','Sun','Earth Position at Launch','Jupiter Position at Arrival','Saturn Position at Arrival')
% Plots the trajectory of probe 2
figure(6)
plot3(Earth_Position_ICRF(:,1),Earth_Position_ICRF(:,2),Earth_Position_ICRF(:,3),'b');grid on; hold on; % Earth orbit in ICRF frame
plot3(Jupiter_Position_ICRF(:,1),Jupiter_Position_ICRF(:,2),Jupiter_Position_ICRF(:,3),'r'); % Jupiter orbit in ICRF frame
plot3(Uranus_Arrival_Position_Vector_ICRF(:,1),Uranus_Arrival_Position_Vector_ICRF(:,2),Uranus_Arrival_Position_Vector_ICRF(:,3),'m'); % Uranus orbit in ICRF frame
plot3(X_Earth_Jupiter_spacecraft_orbit,Y_Earth_Jupiter_spacecraft_orbit,Z_Earth_Jupiter_spacecraft_orbit,'k'); % Earth Jupiter Spacecraft orbit
plot3(X_Jupiter_Uranus_spacecraft_orbit,Y_Jupiter_Uranus_spacecraft_orbit,Z_Jupiter_Uranus_spacecraft_orbit,'c') % Jupiter Uranus Spacecraft orbit
plot3(0,0,0, 'o', 'MarkerFaceColor', 'y') % yellow dot representing sun, not to scale
plot3(Position_Earth_Launch(1),Position_Earth_Launch(2),Position_Earth_Launch(3), 'o', 'MarkerFaceColor', 'b'); % earth position at launch
plot3(Position_Jupiter_Arrival(1),Position_Jupiter_Arrival(2),Position_Jupiter_Arrival(3), 'o', 'MarkerFaceColor', 'r') % jupiter position at arrival
plot3(Uranus_Position_Arrival(1),Uranus_Position_Arrival(2),Uranus_Position_Arrival(3),'o','MarkerFaceColor', 'm')% Uranus position at arrival
title('Probe 2:Earth-Jupiter-Uranus Trajectory in ICRF frame')
xlabel('X [km]');ylabel('Y [km]');zlabel('Z [km]')
legend('Earth Orbit','Jupiter orbit','Uranus Orbit','Earth-Jupiter Spacecraft Coasting Orbit','Jupiter-Uranus Spacecraft Coasting Orbit','Sun','Earth Position at Launch','Jupiter Position at Arrival','Uranus Position at Arrival')
% Plots the combine trajectory of probe 1&2
figure(7)
plot3(Earth_Position_ICRF(:,1),Earth_Position_ICRF(:,2),Earth_Position_ICRF(:,3),'b');grid on; hold on; % Earth orbit in ICRF frame
plot3(Jupiter_Position_ICRF(:,1),Jupiter_Position_ICRF(:,2),Jupiter_Position_ICRF(:,3),'r'); % Jupiter orbit in ICRF frame
plot3(Saturn_Arrival_Position_Vector_ICRF(:,1),Saturn_Arrival_Position_Vector_ICRF(:,2),Saturn_Arrival_Position_Vector_ICRF(:,3),'m'); % Saturn orbit in ICRF frame
plot3(Uranus_Arrival_Position_Vector_ICRF(:,1),Uranus_Arrival_Position_Vector_ICRF(:,2),Uranus_Arrival_Position_Vector_ICRF(:,3),'y'); % Uranus orbit in ICRF frame
plot3(X_Earth_Jupiter_spacecraft_orbit,Y_Earth_Jupiter_spacecraft_orbit,Z_Earth_Jupiter_spacecraft_orbit,'k'); % Earth Jupiter Spacecraft orbit
plot3(X_Jupiter_Saturn_spacecraft_orbit,Y_Jupiter_Saturn_spacecraft_orbit,Z_Jupiter_Saturn_spacecraft_orbit,'g') % Jupiter Saturn Spcecraft orbit
plot3(X_Jupiter_Uranus_spacecraft_orbit,Y_Jupiter_Uranus_spacecraft_orbit,Z_Jupiter_Uranus_spacecraft_orbit,'c') % Jupiter Uranus Spacecraft orbit
plot3(0,0,0, 'o', 'MarkerFaceColor', 'y') % yellow dot representing sun, not to scale
plot3(Position_Earth_Launch(1),Position_Earth_Launch(2),Position_Earth_Launch(3), 'o', 'MarkerFaceColor', 'b'); % earth position at launch
plot3(Position_Jupiter_Arrival(1),Position_Jupiter_Arrival(2),Position_Jupiter_Arrival(3), 'o', 'MarkerFaceColor', 'r') % jupiter position at arrival
plot3(Saturn_Position_Arrival(1),Saturn_Position_Arrival(2),Saturn_Position_Arrival(3),'o', 'MarkerFaceColor', 'm') % saturn position at arrival
plot3(Uranus_Position_Arrival(1),Uranus_Position_Arrival(2),Uranus_Position_Arrival(3),'o','MarkerFaceColor', "#A2142F")% Uranus position at arrival
title('Probe 1:Earth-Jupiter-Saturn Trajectory in ICRF frame and Probe 2:Earth-Jupiter-Uranus Trajectory in ICRF frame')
xlabel('X [km]');ylabel('Y [km]');zlabel('Z [km]')
legend('Earth Orbit','Jupiter orbit','Saturn Orbit','Uranus Orbit','Earth-Jupiter Spacecraft Coasting Orbit','Jupiter-Saturn Spacecraft Coasting Orbit','Jupiter-Uranus Spacecraft Coasting Orbit','Sun','Earth Position at Launch','Jupiter Position at Arrival','Saturn Position at Arrival','Uranus Position at Arrival')
% Outputs of launch and arrival dates and flyby delta-v gain.
Final_Earth_Launch_Date = ['Earth Launch Date: ',string(Achievable_Earth_Launch_Date_string)];disp(Final_Earth_Launch_Date);
Final_Earth_Hyperbolic_Launch_Velocity = ['Hyperbolic Departure Velocity From Earth: ',num2str(Achivable_Earth_Departure_Hyperbolic_Velocity),'km/sec'];disp(Final_Earth_Hyperbolic_Launch_Velocity);
Final_Jupiter_Arrival_Date = ['Jupiter Arrival Date: ',string(Achievable_Jupiter_Arrival_Date_string)];disp(Final_Jupiter_Arrival_Date);
Final_Delta_v_Jupiter_Flyby = ['Delta-v Jupiter Flyby: ',num2str(delta_v_Jupiter_flyby),'km/sec'];disp(Final_Delta_v_Jupiter_Flyby);
Final_Jupiter_Departure_Date = ['Jupiter Departure Date: ',string(Jupiter_Departure_Date)];disp(Final_Jupiter_Departure_Date);
Final_Saturn_Arrival_Date = ['Saturn Arrival Date: ',string(Saturn_Arrival_Date_String)];disp(Final_Saturn_Arrival_Date);
Final_Delta_v_Saturn_Flyby = ['Delta-v Saturn Flyby: ',num2str(delta_v_Saturn_flyby),'km/sec'];disp(Final_Delta_v_Saturn_Flyby);
Final_Saturn_Departure_Date = ['Saturn Departure Date: ',string(Saturn_Arrival_Date_String)];disp(Final_Saturn_Departure_Date);
Final_Uranus_Arrival_Date = ['Uranus Arrival Date: ',string(Uranus_Arrival_Date_String)];disp(Final_Uranus_Arrival_Date);
Final_Delta_v_Uranus_Flyby = ['Delta-v Uranus Flyby: ',num2str(delta_v_Uranus_flyby),'km/sec'];disp(Final_Delta_v_Uranus_Flyby);
Final_Uranus_Departure_Date = ['Uranus Departure Date: ',string(Uranus_Arrival_Date_String)];disp(Final_Uranus_Departure_Date);