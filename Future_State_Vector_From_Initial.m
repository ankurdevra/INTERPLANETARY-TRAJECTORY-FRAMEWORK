function [Future_State_Vectors] = Future_State_Vector_From_Initial(Initial_State_Vector,Final_Time)
% Calculates the final state vectors of a body orbitting around earth
% from the inital position and velocity vector w.r.t ECI frame by employing RK 4 numerical scheme.
% Also plots the motion of body around earth till specified final time in ECI frame.
% REQUIRED INPUTS:
% Initial_State_Vector = row vector of initial state vector of orbitting
% body in km and km/sec
% Final_Time = Final time till which to find the state vectors for in sec
% OUTPUTS:
% Future_State_Vectors = Calculated state vectors at the final time in km
% and km/sec as a column vector
%% Creator:- ANKUR DEVRA 
% Develope Date - 22 March 2022
% Iteration 1 -
%% Starting Data
format long;
y0 = Initial_State_Vector'; % column vector of initial state vectors in km and km/sec
mu_earth=398600.4418; % Gravitational parameter of earth in km^3/s^2
%% RK 4 numerical scheme
h=10;%step size for integration
%initial values
t0 = 0; % sec
tend = Final_Time; % sec, 
Steps = tend/h; % how long to run for loop
% memory pre allocation
X = zeros(length(Steps));
Y = zeros(length(Steps));
Z = zeros(length(Steps));
VX = zeros(length(Steps));
VY = zeros(length(Steps));
VZ = zeros(length(Steps));
% integration from t0 to tfinal using a for loop
for i=1:Steps
    y   = RK4(@Derivative,t0, y0, h);
    y0 = y;
    t0 = t0+h;
    X(i) = y0(1);% x position
    Y(i) = y0(2);% y position
    Z(i) = y0(3);% z position
    VX(i) = y0(4);% x velocity
    VY(i) = y0(5);% y velocity
    VZ(i) = y0(6);% z velocity
end
%% Output
Future_State_Vectors = [X(end) Y(end) Z(end) VX(end) VY(end) VZ(end)]';
%% plotting earth
imData = imread('2_no_clouds_4k.jpg');% loads flat projection of earth
[xS,yS,zS] = sphere(50); % returns the x-, y-, and z- coordinates of a sphere with a radius equal to 1 and 50-by-50 faces.
earth_radius = (6378137.0)*10^-3;  % radius of earth in km
xSE = earth_radius*xS;% scales the sphere in x direction by a factor of earth radius
ySE = earth_radius*yS;% scales the sphere in y direction by a factor of earth radius
zSE = earth_radius*zS;% scales the sphere in z direction by a factor of earth radius
figure(1)
surface(xSE,ySE,zSE);hold on % plots the surface of earth
axis equal;
ch = get(gca,'children'); % get properties of current figure
set(ch,'facecolor','texturemap','cdata',flipud(imData),'edgecolor','none') % asthetic properties
set(gca,'color','black')% asthetic properties
%% plotting orbit
plot3(X,Y,Z,'w','LineWidth',2);
title('orbit in ECI frame');
xlabel('ECI x (km)')
ylabel('ECI y (km)')
zlabel('ECI z (km)')
grid on ;
axis vis3d
hold off
%% RK 4 scheme
    function ya = RK4(Derivative,t0, y0, h)
        % Runga Kutta to solve 2BEOM
        k1 = Derivative(t0    , y0          );
        k2 = Derivative(t0+h/2, y0+(h/2)*k1);
        k3 = Derivative(t0+h/2, y0+(h/2)*k2);
        k4 = Derivative(t0+h  , y0+    h*k3);
        ya = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
%% Differential equation in state space form ( Modified Newton law of universal gravitation )
     function dydt = Derivative(~,y)
        r=sqrt((y(1)^2)+(y(2)^2)+(y(3)^2));% calculated the L2 norm
        dydt=[y(4); y(5); y(6); -(mu_earth/r^3)*y(1); -(mu_earth/r^3)*y(2); -(mu_earth/r^3)*y(3)];% State vector derivative
    end
end