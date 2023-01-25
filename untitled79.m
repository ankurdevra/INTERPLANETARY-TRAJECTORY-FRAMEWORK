%% Starting data
clc;clear
Departure_planet = 'Earth';
Arrival_planet = 'Jupiter';
Year_start = 2036; Month_start=9; Date_start=25;
Year_end = 2039; Month_end=12; Date_end=31;
mu_sun = 1.32712440018*10^(11); % km^3/s^2 sun gravitational parameter
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal parameter
R_earth = 6378.137; % equatorial radius of earth km
% Time frame
JD_start = juliandate([Year_start Month_start Date_start]);
JD_end = juliandate([Year_end Month_end Date_end]);
Time_range_JD = JD_start:1:JD_end;
% Memory preallocation
R_departure_planet = (zeros(3,length(Time_range_JD)))';
V_departure_planet = (zeros(3,length(Time_range_JD)))';
R_arrival_planet = (zeros(3,length(Time_range_JD)))';
V_arrival_planet = (zeros(3,length(Time_range_JD)))';
R_departure_planet_norm = (zeros(1,length(Time_range_JD)))';
V_departure_planet_norm = (zeros(1,length(Time_range_JD)))';
R_arrival_planet_norm = (zeros(1,length(Time_range_JD)))';
V_arrival_planet_norm = (zeros(1,length(Time_range_JD)))';
v_inf = (zeros(1,length(Time_range_JD)))';
delta_v = (zeros(1,length(Time_range_JD)))';
v_inf_arrival = (zeros(1,length(Time_range_JD)))';
for i = 1:length(Time_range_JD)
    [R_departure_planet(i,1:end),V_departure_planet(i,1:end)] = planetEphemeris(Time_range_JD(i),'Sun',Departure_planet,'432t'); % km and km/sec,heliocentric coordinates of departure planet at the
    % time of departure in ICRF coords frame
    R_departure_planet_norm(i) = norm(R_departure_planet(i,1:end)); % km, distance of orbiting departure planet from sun in ICRF coords
    V_departure_planet_norm(i) = norm(V_departure_planet(i,1:end)); % km/sec velocity of orbiting departure planet from sun in ICRF coords

    [R_arrival_planet(i,1:end),V_arrival_planet(i,1:end)] = planetEphemeris(Time_range_JD(i),'Sun',Arrival_planet,'432t');% km and km/sec,heliocentric coordinates of arrival planet at the
    % time of arrival in ICRF coords frame
    R_arrival_planet_norm(i) = norm(R_arrival_planet(i,1:end)); % km, distance of orbiting arrival planet from sun in ICRF coords
    V_arrival_planet_norm(i) = norm(V_arrival_planet(i,1:end)); % km/sec velocity of orbiting arrival planet from sun in ICRF coords


end
% figure(1)
% plot3(Time_range_JD,R_departure_planet_norm,V_departure_planet_norm,'k');hold on;plot3(Time_range_JD,R_arrival_planet_norm,V_arrival_planet_norm,'r');
% grid on;
for i = 1:length(Time_range_JD)
    v_inf(i) = sqrt(mu_sun/(R_departure_planet_norm(i)))*((sqrt((2*R_arrival_planet_norm(i))/(R_departure_planet_norm(i)+R_arrival_planet_norm(i))))-1);% km/sec, hyperbolic excess speed required by spacecraft in ICRF frame
end

rp = 300; % km, radius of initial circular parking orbit 
v_c = sqrt(mu_earth/(R_earth+rp)); % km/sec, speed of spacecraft in its initial parking orbit
for i = 1:length(Time_range_JD)
    delta_v(i) = v_c*(((sqrt(2+(((v_inf(i))/(v_c))^2)))-1)); % km/sec, delta v required to set up the departure hyperbola
end
%
for i = 1:length(Time_range_JD)
    v_inf_arrival(i) = sqrt(mu_sun/(R_arrival_planet_norm(i)))*(1-(sqrt((2*R_departure_planet_norm(i))/(R_departure_planet_norm(i)+R_arrival_planet_norm(i)))));% km/sec, arrival hyperbolic excess speed by spacecraft in ICRF frame
end
%
figure(1)
plot(Time_range_JD,delta_v,'k');hold on;
plot(Time_range_JD,v_inf,'r');
plot(Time_range_JD,v_inf_arrival,'g')
grid on;
xlabel('Julian Dates');ylabel('Departure and arrival Hyperbolic excess speed/Delta v required for departure hyerbola [km/sec]')
legend('Delta v','Departure Hyperbolic excess velocity','Arrival Hyperbolic excess velocity')
