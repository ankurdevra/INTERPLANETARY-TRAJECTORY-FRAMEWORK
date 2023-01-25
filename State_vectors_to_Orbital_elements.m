%% ANKUR DEVRA
% The following code calculates eccentricity, semimajor axis, semilatus
% rectum, time period of orbit, radius of apogee and perigee and specific
% angular momentum,radial velocity,inclination,RAAN,argument of
% perigee,True anomaly,Eccentric anomaly,Mean anomaly,and number of orbit
% in a day.
% from given orbital state vectors. It also plots orbit of satellite for one time peiod
% using RK4
clc;clear;
format long;
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
% Initial state of satellite (x,y,z,vx,vy,vz) in km and km/sec must be
% column matrix
y0 = [-5557351.475/1000,3029589.686/1000,-2608217.221/1000,-4315.974339/1000,-3306.557896/1000,5352.376949/1000]';
%y0=[6513.7755,-2663.0167,1.7061,-1.303,3.94994,7.5397]';
%y0=(10^8).*[1.049985872084391,0.960147117972740,0.416286941925393,-21.514912060056147*10^(-8),19.266421268391959*10^(-8),8.351982541771969*10^(-8)]';
%y0=[-6510.776,-2263.017,14.706,-0.303004,0.949936,7.539667]';
%y0 = [-134883.9813,374965.1607,35659.8733,-0.93188,-0.29528,-0.0084634]';%state vectors of moon

% some extra state vectors.
%{ 
%y0 = [4095.0191,1007.1462,5328.0591,-1.7370262,7.4617554,-0.0754351]';% state vector of ISS
%y0 = [7815820.187/1000,-14555684.988/1000,-20700198.639/1000,3717.547/1000,648.838/1000,922.737/1000]';
%} 
r_vec = y0(1:3);% distance vector in km 
r_mag = norm(r_vec);% initial distance in km
v_vec = y0(4:6);% speed vector in km/sec
v_mag = norm(v_vec);% initila velocity in km/sec
h_vec = cross(r_vec,v_vec);% specific angular momentum vector km^2/sec
h_mag = norm(h_vec); % specific angular momentum km^2/sec
radial_velocity = (dot(r_vec,v_vec))/r_mag; %km/sec (v_r) if v_r>0 trajectory away from preigee if v_r<0 towards perigee
e_vec = (1/mu_earth)*(((v_mag^2)-(mu_earth/r_mag))*r_vec - r_mag*radial_velocity*v_vec); % eccentricity vector
e_mag = norm(e_vec); % eccentricity
T = (((2*pi)/(mu_earth^2))*((h_mag/(sqrt(1-e_mag^2))))^3)/3600; % Time period of orbit in hrs
rp = ((h_mag^2)/(mu_earth))*(1/(1+e_mag*cosd(0))); % radius of perigee in km
ra = ((h_mag^2)/(mu_earth))*(1/(1+e_mag*cosd(180))); % radius of apogee in km
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

semimajor_axis = ['The semi major axis is ',num2str(a),' km'];disp(semimajor_axis);
eccentricity = ['The eccentricity is ',num2str(e_mag),''];disp(eccentricity);
semilatus_rectum = ['The semilatus_rectum is ',num2str(l),' km'];disp(semilatus_rectum);
Time_period = ['The Time Period is ',num2str(T),' hrs'];disp(Time_period);
radius_perigee = ['The radius of perigee is ',num2str(rp),' km'];disp(radius_perigee);
radius_apogee = ['The radius of apogee is ',num2str(ra),' km'];disp(radius_apogee);
if i<90 % prograde orbit
    inclination = ['The inclination of orbit is ',num2str(i),' degrees, prograde orbit'];disp(inclination);
elseif i==90 % polar orbit
    inclination = ['The inclination of orbit is ',num2str(i),' degrees, polar orbit'];disp(inclination);
elseif i>90 % retrograde orbit
    inclination = ['The inclination of orbit is ',num2str(i),' degrees, retrograde orbit'];disp(inclination);
end
RAAN = ['The Right Ascension of Ascending Node of orbit is ',num2str(omega),' degrees'];disp(RAAN);
Argument_of_Perigee = ['The argument of perigee of orbit is ',num2str(small_omega),' degrees'];disp(Argument_of_Perigee);
True_anomaly = ['The True anomaly orbit is ',num2str(theta),' degrees'];disp(True_anomaly);
Eccentric_anomaly = ['The Eccentric anomaly orbit is ',num2str(E),' degrees'];disp(Eccentric_anomaly);
Mean_anomaly = ['The Mean anomaly orbit is ',num2str(M),' degrees'];disp(Mean_anomaly);
Orbit_per_day = ['Number of Orbits per day are ',num2str(Orbit_number),' orbits'];disp(Orbit_per_day);

h=10;%step size for integration 
%initial values
t0 = 0; % sec
tend = T*3600*1.5;% sec, iterating till more than the end of one time period
Steps = tend/h; % how long to run for loop
% memory pre allocation
X = zeros(length(Steps));
Y = zeros(length(Steps));
Z = zeros(length(Steps));
% integration from t0 to tfinal using a for loop 
for i=1:Steps
    y   = RK4(@deriv,t0, y0, h);
    y0 = y;
    t0 = t0+h;
    X(i) = y0(1);
    Y(i) = y0(2);
    Z(i) = y0(3);
end
%plotting earth 
imData = imread('2_no_clouds_4k.jpg');
[xS,yS,zS] = sphere(50);
earth_radius = (6378137.0)*10^-3;  % kilo meters
xSE = earth_radius*xS;
ySE = earth_radius*yS;
zSE = earth_radius*zS;
figure(1)
surface(xSE,ySE,zSE);hold on
axis equal;
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',flipud(imData),'edgecolor','black')
set(gca,'color','black')

% plotting orbit
plot3(X,Y,Z,'w','LineWidth',2);
title('orbit in ECI frame');
xlabel('ECI x (km)')
ylabel('ECI y (km)')
zlabel('ECI z (km)')
grid on ;
axis vis3d
hold off
%view([42.8864468,-78.8783689,45])
% to rotate plot of orbit in 3D
%{
%view([42.8864468,-78.8783689,45]) % the x and y values are latitude and longitude of BUFFALO, NY
%line
% for AZ = 0:0.3:360
%   view(AZ, 23.5); % rotation about earth axial tilt
%   pause(0.01);
% end
% for EL = -90:2:90
%   view(45, EL);
%   pause(0.1);
% end
%}
function ya = RK4(deriv,t0, y0, h)
% Runga Kutta to solve 2BEOM
k1 = deriv(t0    , y0          );
k2 = deriv(t0+h/2, y0+(h/2)*k1);
k3 = deriv(t0+h/2, y0+(h/2)*k2);
k4 = deriv(t0+h  , y0+    h*k3);
ya = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end
function dydt = deriv(~,y)
% State vector derivative
    mu_earth=398600.4418; %km^3/s^2
    r=sqrt((y(1)^2)+(y(2)^2)+(y(3)^2));% calculated the L2 norm
    dydt=[y(4); y(5); y(6); -(mu_earth/r^3)*y(1); -(mu_earth/r^3)*y(2); -(mu_earth/r^3)*y(3)];
end