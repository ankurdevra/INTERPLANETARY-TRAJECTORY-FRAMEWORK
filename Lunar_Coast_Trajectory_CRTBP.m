function [Spacecraft_Position_and_Velocity,Jacobi_Constant,Initial_Specific_Angular_Momentum] = Lunar_Coast_Trajectory_CRTBP(Initial_Altitude_of_Spacecraft_above_Earth,Azimuth,Initial_Flight_Path_Angle,Relative_Burnout_Velocity,Time_days)
% Calculates the position and velocity of a spacecraft coasting under the
% influence of earth moon gravity system. Plots the translunar coasting of the
% spacecraft with respect to the barycenter reference frame.
% Also gives the jacobi constant and initial specific angular momentum of
% spacecraft.
% Barycenter reference frame is non inertial (comoving frame)
% Accounted for the inertial angular velocity due to non inertial frame.
% Barycenter coordinates (0,0,0)
% Assumes planar coasting no z position/velocity component.
% REQUIRED INPUTS:
% Initial_Altitude_of_Spacecraft_above_Earth = Initial Radial altitude
% above earth km
% Azimuth = planar azimuth w.r.t earth moon line degrees
% Initial_Flight_Path_Angle = Initial flight path angle of spacecraft in
% deg
% Relative_Burnout_Velocity = Burnout velocity of spacecraft km/sec wrt
% barycenter
% Time_days = How long in days to propagate the orbit for; days
% OUTPUTS:
% Spacecraft_Position_and_Velocity = final state vector of spacecraft wrt
% barycenter in km and km/sec
% Jacobi_Constant = Jacobi constant of spacecraft, inidication of total
% energy relative to barycenter frame; km^2/sec^2
% Initial_Specific_Angular_Momentum = Initial specific angular momentum of
% spacecraft km^2/sec
% Sample starting point:[Spacecraft_Position_and_Velocity,Jacobi_Constant,Initial_Specific_Angular_Momentum] = Lunar_Coast_Trajectory_CRTBP(200,-90,20,10.915,4)
%% Creator:- ANKUR DEVRA 
% Develope Date - 23 March 2022
% Iteration 1 -
%% Starting data
format long;
mass_earth = 5.9722 * 10^24; % kg mass of earth 
mass_moon = 0.07346 * 10^24; % kg mass of moon
R_earth = 6378; % radius of earth km
R_moon = 1738.1; % raius of moon km
G = 6.67430*10^(-20); % universal gravitational constant km^3/kg-sec^2
r12 = 384399; % km distance between earth and moon
m1 = mass_earth; % kg
m2 = mass_moon; % kg
pi_1 = m1/(m1+m2); % mass ratio (dimensionless)
pi_2 = m2/(m1+m2); % mass ratio ( dimensionless)
x1 = -pi_2*r12; % km x coordinate of earth (relative to barycenter)
x2 = pi_1*r12; % km x coordinate of moon (relative to barycenter)
mu1 = G*m1;% Gravitational parameter of earth in km^3/s^2
mu2 = G*m2;% Gravitational parameter of moon in km^3/s^2
mu_system = G*(m1+m2); % Gravitational parameter of earth moon system
d = Initial_Altitude_of_Spacecraft_above_Earth; % radial distance of spacecraft above earth km
phi = Azimuth; % initial azimuth of spacecraft deg
v_bo = Relative_Burnout_Velocity; % burnout velocity of spacecraft km/sec
Final_Time = Time_days*86400; % days converted to sec; sec
gamma = Initial_Flight_Path_Angle; % initial flight path angle of spacecraft deg
x = x1 + (d+R_earth)*cosd(phi); % km form trig; x coordinate of spacecraft from barycenter reference frame
y = (d+R_earth)*sind(phi); % km form trig; y coordinate of spacecraft from barycenter reference frame
z = 0; % km z coordinate of spacecraft from barycenter reference frame
%vx = v_bo*(sind(gamma-phi)); %km/sec x component of spacecraft velocity from barycenter reference frame 
%vy = v_bo*(sind(gamma+phi)); %km/sec y component of spacecraft velocity from barycenter reference frame 
vx = v_bo*(sind(gamma)*cosd(phi) - cosd(gamma)*sind(phi));
vy = v_bo*(sind(gamma)*sind(phi) + cosd(gamma)*cosd(phi));
vz = 0; %km/sec z component of spacecraft velocity from barycenter reference frame 
Omega = sqrt(((mu_system))/((r12)^3)); % inertial angular velocity due to non inertial frame rad/sec
y0 = [x,y,z,vx,vy,vz]';% Initial state vector w.r.t barycenter reference frame km and km/sec
r1=sqrt((x + (pi_2*r12))^2 + y^2 + z^2);% km initial distance of space craft from earth
r2=sqrt((x - (pi_1*r12))^2 + y^2 + y^2);% km initial distance of space craft from moon
C = ((v_bo)^2)/2 - (1/2)*((Omega)^2)*((x)^2 + (y)^2) - ( mu1/r1) -( mu2/r2);% calculated jacobi constant
h = (d+R_earth)*v_bo*cosd(phi); % km^2/sec initial angular momentum of spacecraft
%% RK 4 numerical scheme
h_=10;%step size for integration
%initial values
t0 = 0; % sec
tend = Final_Time; % sec, 
Steps = tend/h_; % how long to run for loop
% memory pre allocation
x = zeros(length(Steps));
y = zeros(length(Steps));
z = zeros(length(Steps));
vx = zeros(length(Steps));
vy = zeros(length(Steps));
vz = zeros(length(Steps));
%integration from t0 to tfinal using a for loop
for i=1:Steps
    y_   = RK4(@Derivative,t0, y0, h_);
    y0 = y_;
    t0 = t0+h_;
    x(i) = y0(1);% x position
    y(i) = y0(2);% y position
    z(i) = y0(3);% z position
    vx(i) = y0(4);% x velocity
    vy(i) = y0(5);% y velocity
    vz(i) = y0(6);% z velocity
end
%% Output
Spacecraft_Position_and_Velocity = [x(end) y(end) z(end) vx(end) vy(end) vz(end)]';
Jacobi_Constant = C;
Initial_Specific_Angular_Momentum = h;
%% plotting earth
figure(1)
imData = imread('2_no_clouds_4k.jpg');% loads flat projection of earth
[xS,yS,zS] = sphere(50); % returns the x-, y-, and z- coordinates of a sphere with a radius equal to 1 and 50-by-50 faces.
earth_radius = R_earth;  % radius of earth in km
xSE = earth_radius*xS;% scales the sphere in x direction by a factor of earth radius
ySE = earth_radius*yS;% scales the sphere in y direction by a factor of earth radius
zSE = earth_radius*zS;% scales the sphere in z direction by a factor of earth radius
surface(xSE+x1,ySE,zSE);hold on % plots the surface of earth
axis vis3d
ch = get(gca,'children'); % get properties of current figure
set(ch,'facecolor','texturemap','cdata',flipud(imData),'edgecolor','none') % asthetic properties
set(gca,'color','white')% asthetic properties
%% plotting moon
imDatam = imread('1_moon_flat_projection.jpg');% loads flat projection of moon
[xSm,ySm,zSm] = sphere(50); % returns the x-, y-, and z- coordinates of a sphere with a radius equal to 1 and 50-by-50 faces.
moon_radius = R_moon;  % radius of earth in km
xSEm = moon_radius*xSm;% scales the sphere in x direction by a factor of earth radius
ySEm = moon_radius*ySm;% scales the sphere in y direction by a factor of earth radius
zSEm = moon_radius*zSm;% scales the sphere in z direction by a factor of earth radius
surface(xSEm+x2,ySEm,zSEm);% plots the surface of moon
axis vis3d
ch = get(gca,'children'); % get properties of current figure
set(ch,'facecolor','texturemap','cdata',flipud(imDatam),'edgecolor','none') % asthetic properties
set(gca,'color','white')% asthetic properties
%% plotting translunar orbit
plot3(x,y,z,'k','LineWidth',1);
title('Translunar Coasting of Spacecraft');
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
grid on ;
axis([(min(x)-20000) (max(x)+20000) (min(y)-10000) (max(y)+10000) -50000 50000])
axis vis3d
hold off
%% Animates the combined motion of both the masses
% figure(2)
% imData = imread('2_no_clouds_4k.jpg');% loads flat projection of earth
% [xS,yS,zS] = sphere(50); % returns the x-, y-, and z- coordinates of a sphere with a radius equal to 1 and 50-by-50 faces.
% earth_radius = R_earth;  % radius of earth in km
% xSE = earth_radius*xS;% scales the sphere in x direction by a factor of earth radius
% ySE = earth_radius*yS;% scales the sphere in y direction by a factor of earth radius
% zSE = earth_radius*zS;% scales the sphere in z direction by a factor of earth radius
% surface(xSE+x1,ySE,zSE);hold on % plots the surface of earth
% axis vis3d
% ch = get(gca,'children'); % get properties of current figure
% set(ch,'facecolor','texturemap','cdata',flipud(imData),'edgecolor','none') % asthetic properties
% set(gca,'color','white')% asthetic properties
% plot3(0,0,0,'ko','MarkerFaceColor','k');hold on % inertial origin
% title('Motion of Bodies relative to the Inertial Reference Frame Animation');
% xlabel('Inertial X (km)');
% ylabel('Inertial Y (km)');
% zlabel('Inertial Z (km)');
% axis equal
% for i = 2:Steps
%     plot3(x(i),y(i),z(i),'r');hold on;grid on % body 1
%     plot3([x(i-1) x(i)],[y(i-1) y(i)],[z(i-1) z(i)],'r');
%     pause(0.001)
% end
%% RK 4 scheme
    function ya = RK4(Derivative,t0, y0, h)
        % Runga Kutta to solve CRTBP
        k1 = Derivative(t0    , y0          );
        k2 = Derivative(t0+h/2, y0+(h/2)*k1);
        k3 = Derivative(t0+h/2, y0+(h/2)*k2);
        k4 = Derivative(t0+h  , y0+    h*k3);
        ya = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
%% Differential equation in state space form ( Modified Newton law of universal gravitation )
     function dydt = Derivative(~,y)
        r1_=sqrt((y(1) + (pi_2*r12))^2 + y(2)^2 + y(3)^2);% km  distance of space craft from earth
        r2_=sqrt((y(1) - (pi_1*r12))^2 + y(2)^2 + y(3)^2);% km  distance of space craft from moon
        % dydt contain velocity and acceleration of spacecraft 
        dydt=[y(4); y(5); y(6); 2*Omega*y(5) + y(1)*Omega^2 - (mu1/r1_^3)*(y(1)+(pi_2*r12)) - (mu2/r2_^3)*(y(1)-(pi_2*r12));-2*Omega*y(4) + y(2)*Omega^2 - (mu1/r1_^3)*(y(2)) - (mu2/r2_^3)*(y(2)); -(mu1/r1_^3)*(y(3)) - (mu2/r2_^3)*(y(3))];% State vector derivative
    end
end