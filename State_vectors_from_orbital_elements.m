%% ANKUR DEVRA
% The following code converts orbital state elements to orbital state vectors in
% perfocal frame, geocentric equatorial frame and ECI frame. (geocentric equatorial and ECI are same frames only)
%clc;clear;format long;

% UNSTABLE CODE, NEEDS REFINEMENT, INSTEAD USE
% "Geocentric_equi_and_perifocal_state_vector_from_orbital_element"
% FUNCTION

function [Perifocal_coordinates,Geocentric_equatorial_coordinates] = State_vectors_from_orbital_elements(Semimajor_axis,Eccentricity,Inclination,Argument_of_perigee,Mean_anomaly,RAAN)
%% data for moon
a = Semimajor_axis; % semi major axis km
e = Eccentricity; % eccentricity of moon
i = Inclination; % degrees inclination of moon 
small_omega = Argument_of_perigee; % degress argument of perigee of moon
M = Mean_anomaly; % degrees mean anomaly of moon
Omega = RAAN; % degress RAAN
%%
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
h = sqrt(a*mu_earth*(1-(e^2))); %km^2/sec angular momentum magnitude
M_rad = M*(pi/180);% radian mean anomaly
T = (((2*pi)/sqrt(mu_earth))*((a)^(3/2)))/3600; % Time period of orbit in hrs
% used newton iteration to solve for eccentric and true anomaly
if M_rad > pi
    E = M_rad - e/2; % initial guess
elseif M_rad == pi
    E = M;
elseif M_rad < pi
    E = M_rad + e/2; % initial guess
end

matrix=[]; % initailize empty matrix
while abs(((E-e*sin(E)-M_rad)/(1-e*cos(E)))) > 10^(-8) % specifying error tolerance
    E = E - ((E-e*sin(E)-M_rad)/(1-e*cos(E))); % ECCENTRIC ANOMALY USING NEWTON iteration METHOD
    theta = 2*atan((sqrt((1+e)/(1-e)))*tan(E/2));%+360; % true anomaly for each value of eccentric anomaly
    if theta<0
        matrix = [matrix;[E,theta+2*pi]];% stores the result in a matrix after each iteration
    else
        matrix = [matrix;[E,theta]];% stores the result in a matrix after each iteration
    end
end

eccentric_anomaly = matrix(:,1);
true_anomaly = matrix(:,2);
eccentric_anomaly_deg = (matrix(end,1));%*(180/pi);% eccentric anomaly in degrees % keep this in radian as all subsequent calculation of angle works for radian
true_anomaly_deg = (matrix(end,2));%*(180/pi);% true anomaly in degrees % keep this in radian as all subsequent calculation of angle works for radian
theta = true_anomaly_deg; % true anomaly in degrees
perifocal_position_vec = [cos(theta);sin(theta);0];%transform matrix
perifocal_velocity_vec = [-sin(theta);e+cos(theta);0];%transform matrix
r_x_bar = (((h^2)/(mu_earth))*(1/(1+e*cos(theta)))).*perifocal_position_vec; % km position vector in perifocal coordinates
v_x_bar = (mu_earth/h).*perifocal_velocity_vec; % km/sec velocity vector in perifocal coordinates

Perifocal_coordinates = [r_x_bar,v_x_bar]; % km and km/sec perifocal state vectors

% perifocal to geocentric equatorial transformation using classical euler angle sequence transformation:-
% [Q]_X_x_bar = [R3(arg. of perigee)] [R1(inclination)] [R3(RAAN)] DCM from geocentric equatorial to perifocal
R3_small_omega = [cosd(small_omega) sind(small_omega) 0;-sind(small_omega) cosd(small_omega) 0; 0 0 1]; % degrees arg of perigee rotation
R1_i = [1 0 0;0 cosd(i) sind(i); 0 -sind(i) cosd(i)];% degrees inclination rotation
R3_Omega = [cosd(Omega) sind(Omega) 0;-sind(Omega) cosd(Omega) 0;0 0 1];% degrees RAAN rotation
Q_X_x_bar = (R3_small_omega)*(R1_i)*(R3_Omega);%DCM for transformation from geocentric to perifocal frame degrees
Q_x_bar_X = Q_X_x_bar'; % DCM from perifocal to geocentric degrees
r_X = Q_x_bar_X*r_x_bar; % km position vector in geocentric equatorial frame
v_X = Q_x_bar_X*v_x_bar; % km/sec velocity vector in geocentric equatorial frame

% Perifocal to ECI transformation using classical euler angle sequence:-
% [Q]_perifocal_to_ECI = [R3(RAAN)] [R1(inclination)] [R3(arg. of perigee)] DCM from  perifocal to ECI frame
R3_Omega = [cosd(Omega) -sind(Omega) 0;sind(Omega) cosd(Omega) 0; 0 0 1]; % degrees RAAN rotation
R1_i = [1 0 0;0 cosd(i) -sind(i); 0 sind(i) cosd(i)];% degrees inclination rotation
R3_small_omega = [cosd(small_omega) -sind(small_omega) 0;sind(small_omega) cosd(small_omega) 0; 0 0 1]; % degrees arg of perigee rotation
Q_perifocal_to_ECI = (R3_Omega)*(R1_i)*(R3_small_omega);%DCM for transformation from perifocal to ECI frame degrees
r_ECI = Q_perifocal_to_ECI*r_x_bar; % km position vector in ECI frame
v_ECI = Q_perifocal_to_ECI*v_x_bar; % km/sec velocity vector in ECI frame

%NOTE:- Geocentric equatorial and ECI frames are the same.

Perifocal_x = r_x_bar(1); % x position coordinate in perifocal frames
Perifocal_y = r_x_bar(2); % y position coordinate in perifocal frames
Perifocal_z = r_x_bar(3); % z position coordinate in perifocal frames
Perifocal_x_dot = v_x_bar(1); % x velocity coordinate in perifocal frames
Perifocal_y_dot = v_x_bar(2); % y velocity coordinate in perifocal frames
Perifocal_z_dot = v_x_bar(3); % z velocity coordinate in perifocal frames
P = [Perifocal_x,Perifocal_y,Perifocal_z];
r=norm(P);
ECI_x = r_ECI(1);% x position coordinate in ECI frames
ECI_y = r_ECI(2);% y position coordinate in ECI frames
ECI_z = r_ECI(3);% z position coordinate in ECI frames
ECI_x_dot = v_ECI(1); % x velocity coordinate in ECI frame
ECI_y_dot = v_ECI(2); % y velocity coordinate in ECI frame
ECI_z_dot = v_ECI(3); % z velocity coordinate in ECI frame

ECI_state_vector = [ECI_x,ECI_y,ECI_z,ECI_x_dot,ECI_y_dot,ECI_z_dot]'; % state vectors in ECI frame

Geocentric_equatorial_coordinates = ECI_state_vector; % km and km/sec state vector in geocentric equatorial frame

% Perfocal_coordinates_position = ['The perifocal position coordinates are in km: ',num2str(Perifocal_x),'p + (',num2str(Perifocal_y),')q + ',num2str(Perifocal_z),'w'];disp(Perfocal_coordinates_position)
% Perfocal_coordinates_velocity = ['The perifocal velocity coordinates are in km/sec: ',num2str(Perifocal_x_dot),'p + (',num2str(Perifocal_y_dot),')q + ',num2str(Perifocal_z_dot),'w'];disp(Perfocal_coordinates_velocity);disp(' ');
% 
% ECI_coordinates_position = ['The ECI position coordinates are in km: ',num2str(ECI_x),'i + (',num2str(ECI_y),')j + ',num2str(ECI_z),'k'];disp(ECI_coordinates_position)
% ECI_coordinates_velocity = ['The ECI velocity coordinates are in km/sec: ',num2str(ECI_x_dot),'i + (',num2str(ECI_y_dot),')j + ',num2str(ECI_z_dot),'k'];disp(ECI_coordinates_velocity);disp(' ');
% 
% Time_period = ['The Time Period is ',num2str(T),' hrs'];disp(Time_period);disp(' ');
% 
% y0 = ['y0 = [',num2str(ECI_x),',',num2str(ECI_y),',',num2str(ECI_z),',',num2str(ECI_x_dot),',',num2str(ECI_y_dot),',',num2str(ECI_z_dot),']'';'];disp(y0) % state vectors in ECI frame km and km/sec

end

