function [Angular_momentum,Inclination,Eccentricity,RAAN,Argument_of_Perigee,True_anomaly,Radial_velocity,Time_period_of_orbit_hrs,Radius_perigee,Radius_apogee,Semimajor_axis,Semilatus_rectum,Eccentric_anomaly,Mean_anomaly,Orbits_in_a_day] = Heliocentric_Orbital_elements_from_State_vectors(state_vector)
% ANKUR DEVRA MAE 460
% The following code calculates heliocentric orbit eccentricity, semimajor axis, semilatus
% rectum, time period of orbit, radius of apogee and perigee and specific
% angular momentum,radial velocity,inclination,RAAN,argument of
% perigee,True anomaly,Eccentric anomaly,Mean anomaly,and number of orbit
% in a day from given orbital state vectors.
% SINCE HELIOCENTRIC FRAME IS USED, ORBITAL DATA WILL BE WRT HELIOCENTRIC
% ELLIPTIC FRAME
% REQUIRED INPUTS:
% state_vector = [1x6] km and km/sec heliocentric state vector of an body orbiting sun
% OUTPUTs:
% Angular_momentum = specific angular momentum km^2/sec
% Inclination = inclination of orbit is degrees, if i>90 retrograde orbit
% Eccentricity =  eccentricity
% RAAN =  Right Ascension of Ascending Node, degrees
% Argument_of_Perigee =  argument of perigee degrees
% True_anomaly =  true anomaly degrees
% Radial_velocity =  km/sec (v_r) if v_r>0 trajectory away from preigee if v_r<0 towards perigee
% Time_period_of_orbit_hrs = Time period of orbit in hrs
% Radius_perigee =  radius of perigee in km
% Radius_apogee =  radius of apogee in km
% Semimajor_axis =  semi-major axis in km
% Semilatus_rectum =  semi-latus rectum in km
% Eccentric_anomaly =  Eccentric anomaly in degrees
% Mean_anomaly =  Mean anomaly in degrees
% Orbits_in_a_day =  number of orbits in a day 
%% Creator:- ANKUR DEVRA 
% Develope Date - 2 July 2022
% Iteration 1 -
%% Starting data
mu_sun = 1.32712440018*10^(11);%132712440041.93938;%1.32712440018*10^(11); % km^3/s^2 sun gravitational parameter
%% INPUTS
y0=state_vector'; % initial state vector km and km/sec
%% CALCULATIONS
r_vec = y0(1:3);% distance vector in km
r_mag = norm(r_vec);% initial distance in km
v_vec = y0(4:6);% speed vector in km/sec
v_mag = norm(v_vec);% initial velocity in km/sec
h_vec = cross(r_vec,v_vec);% specific angular momentum vector km^2/sec
h_mag = norm(h_vec); % specific angular momentum km^2/sec
radial_velocity = (dot(r_vec,v_vec))/r_mag; %km/sec (v_r) if v_r>0 trajectory away from preigee if v_r<0 towards perigee
e_vec = (1/mu_sun)*(((v_mag^2)-(mu_sun/r_mag))*r_vec - r_mag*radial_velocity*v_vec); % eccentricity vector
e_mag = norm(e_vec); % eccentricity
T = (((2*pi)/(mu_sun^2))*((h_mag/(sqrt(1-e_mag^2))))^3)/3600; % Time period of orbit in hrs
rp = ((h_mag^2)/(mu_sun))*(1/(1+e_mag*cosd(0))); % radius of perigee in km
ra = ((h_mag^2)/(mu_sun))*(1/(1+e_mag*cosd(180))); % radius of apogee in km
a = (ra+rp)/2; % semi-major axis in km
l = a*(1-(e_mag)^2); % semi-latus rectum in km
h_z = h_vec(3);% z component of specific angular momentum km^2/sec
i = acosd(h_z/h_mag); % inclination of orbit is degrees, if i>90 retrograde orbit
K_cap = [0 0 1]'; % k component vector
Node_vec = cross(K_cap,h_vec);% Node line vector km^2/sec
Node_mag = norm(Node_vec);% Node line magnitude km^2/sec
N_x = Node_vec(1); % x component of Node vector km^2/sec
N_y = Node_vec(2); % y component of Node vector km^2/sec
if N_y >= 0
    omega = acosd(N_x/Node_mag); % Right Ascension of Ascending Node, degrees
elseif N_y < 0
    omega = 360 - acosd(N_x/Node_mag); % Right Ascension of Ascending Node, degrees
end
e_z = e_vec(3); % z component of eccentricity vector
if e_z >= 0 
    small_omega = acosd((dot(Node_vec,e_vec))/(Node_mag*e_mag)); % argument of perigee degrees
elseif e_z < 0
    small_omega = 360 - acosd((dot(Node_vec,e_vec))/(Node_mag*e_mag)); % argument of perigee degrees
end
if radial_velocity >= 0
    theta = acosd((dot(e_vec,r_vec))/(e_mag*r_mag)); %true anomaly degrees
elseif radial_velocity < 0
    theta = 360 - acosd((dot(e_vec,r_vec))/(e_mag*r_mag)); %true anomaly degrees
end
E = (2*atan(((sqrt((1-e_mag)/(1+e_mag)))*tan((theta*(pi/180))/2))))*(180/pi); % Eccentric anomaly in degrees
% converts ecentricity anomlay to positive if necessary
if E<0
    E = E+360;
else
end
M = ((E*(pi/180)) - (e_mag*sin((E*(pi/180)))))*(180/pi); % Mean anomaly in degrees
Orbit_number = (24/T); % number of orbits in a day 
%% OUTPUT:
Angular_momentum = h_mag; % specific angular momentum km^2/sec
Inclination = i; % inclination of orbit is degrees, if i>90 retrograde orbit
Eccentricity = e_mag; % eccentricity
RAAN = omega; % Right Ascension of Ascending Node, degrees
Argument_of_Perigee = small_omega; % argument of perigee degrees
True_anomaly = theta; % true anomaly degrees
Radial_velocity = radial_velocity; % km/sec (v_r) if v_r>0 trajectory away from preigee if v_r<0 towards perigee
Time_period_of_orbit_hrs = T; % Time period of orbit in hrs
Radius_perigee = rp; % radius of perigee in km
Radius_apogee = ra; % radius of apogee in km
Semimajor_axis = a; % semi-major axis in km
Semilatus_rectum = l; % semi-latus rectum in km
Eccentric_anomaly = E; % Eccentric anomaly in degrees
Mean_anomaly = M; % Mean anomaly in degrees
Orbits_in_a_day = Orbit_number; % number of orbits in a day 
end