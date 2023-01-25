function [R1,V1,R2,V2] = Two_Body_Motion_Calculator(Body_1_Mass,Body_2_Mass,Initial_Time,Final_Time,Initial_Position_Body_1,Initial_Position_Body_2,Initial_Velocity_Body_1,Initial_Velocity_Body_2)
% Calculates the final state vectors of two bodies from there inital poition
% and velocity vectors provided there mass and time range is known by employing RK 4 numerical scheme.
% Also plots the motion of both the particle for the
% specified time period and animates the motion.
% REQUIRED INPUTS:
% Body_1_Mass = mass of body 1 in kg
% Body_2_Mass = mass of body 2 in kg
% Initial_Time = Initial time in sec
% Final_Time = Final time in sec
% Initial_Position_Body_1  = Initial position row vector of cartesian position
% for body 1 in km
% Initial_Position_Body_2  = Initial position row vector of cartesian position
% for body 2 in km
% Initial_Velocity_Body_1 = Initial velocity row vector of cartesian
% velocity of body 1 in km/sec
% Initial_Velocity_Body_2 = Initial velocity row vector of cartesian
% velocity of body 2 in km/sec
% OUTPUTS:
% R1 = Final position vector of body 1
% V1 = Final Velocity vector of body 1
% R2 = Final position vector of body 2
% V2 = Final Velocity vector of body 2
% sample starting point:[R1,V1,R2,V2] = Two_Body_Motion_Calculator(10^26,10^26,0,500,[0 0 0],[3000 0 0],[12 20 30],[0 40 0])
% sample starting point:[R1,V1,R2,V2] = Two_Body_Motion_Calculator(10^24,10^26,0,500,[0 0 0],[3000 0 0],[12 20 30],[0 40 0])
% sample starting point:[R1,V1,R2,V2] = Two_Body_Motion_Calculator(10^26,10.3402^26,0,500,[0 0 0],[3000 0 0],[12 20 30],[0 40 0])
% sample starting point:[R1,V1,R2,V2] = Two_Body_Motion_Calculator(10^26,10.056^27,0,500,[0 0 0],[3000 0 0],[10 40 90],[0 40 0])
%% Creator:- ANKUR DEVRA 
% Develope Date - 22 March 2022
% Iteration 1 - 
%% Starting data
format long;
G = 6.67430*10^(-20); % universal gravitational constant km^3/kg-sec^2
R1_0 = Initial_Position_Body_1; % initial position vector of body 1; row vector in km
R2_0 = Initial_Position_Body_2; % initial position vector of body 2; row vector in km
V1_0 = Initial_Velocity_Body_1; % initial velocity vector of body 1; row vector in km/sec
V2_0 = Initial_Velocity_Body_2; % initial velocity vector of body 2; row vector in km/sec
m1 = Body_1_Mass; % mass of body 1 in kg
m2 = Body_2_Mass; % mass of body 2 in kg
Total_mass = m1+m2; % total mass of bodies in kg
y0 = [R1_0,R2_0,V1_0,V2_0]'; % inital state vector of system containing positions and velocities of both the bodies; column vector km and km/sec
%% RK 4 numerical scheme
h=1;%step size for integration 
%initial values
t0 = Initial_Time; % sec starting time of RK4
tend = Final_Time;% sec ending time of RK4 
Steps = (tend/h); % how long to run for loop
% memory pre allocation
X1 = zeros(length(Steps));
Y1 = zeros(length(Steps));
Z1 = zeros(length(Steps));
X2 = zeros(length(Steps));
Y2 = zeros(length(Steps));
Z2 = zeros(length(Steps));
VX1 = zeros(length(Steps));
VY1 = zeros(length(Steps));
VZ1 = zeros(length(Steps));
VX2 = zeros(length(Steps));
VY2 = zeros(length(Steps));
VZ2 = zeros(length(Steps));
% integration from initial time to final time  using a for loop 
for i=1:Steps
    y   = RK4(@Derivative, t0, y0, h);
    y0 = y;
    t0 = t0+h;
    X1(i) = y0(1); % x position body 1
    Y1(i) = y0(2); % y position body 1
    Z1(i) = y0(3); % z position body 1
    X2(i) = y0(4); % x position body 2
    Y2(i) = y0(5); % y position body 2
    Z2(i) = y0(6); % z position body 2
    VX1(i) = y0(7); % x velocity body 1
    VY1(i) = y0(8); % y velocity body 1
    VZ1(i) = y0(9); % z velocity body 1
    VX2(i) = y0(10); % x velocity body 2
    VY2(i) = y0(11); % y velocity body 2
    VZ2(i) = y0(12); % z velocity body 2  
end
%% OUTPUTS
R1 = [X1(end),Y1(end),Z1(end)]'; % final position vector of body 1
R2 = [X2(end),Y2(end),Z2(end)]'; % final position vector of body 2
V1 = [VX1(end),VY1(end),VZ1(end)]'; % final velocity vector of body 1
V2 = [VX2(end),VY2(end),VZ2(end)]'; % final velocity vector of body 2
%% Center of mass calculation
% initializing empty matrices
X_CM = [];
Y_CM = [];
Z_CM = [];
%X_CM = zeros(length(Steps));
for i = 1:(Steps)
    X_CM = [X_CM;(((m1*X1(i)) + (m2*X2(i)))/(Total_mass))];
    Y_CM = [Y_CM;(((m1*Y1(i)) + (m2*Y2(i)))/(Total_mass))];
    Z_CM = [Z_CM;(((m1*Z1(i)) + (m2*Z2(i)))/(Total_mass))];
end
%% Creates random figure # so that there is no over riding of multiple figures
%figure_number = randi([1 1000],1,3);
%% Plots the combined motion of both the mass for the given time frame
figure (1);
plot3(X1,Y1,Z1,'r');hold on; % body 1
plot3(X2,Y2,Z2,'g'); % body 2
plot3(X_CM,Y_CM,Z_CM,'b'); % plots center of mass
plot3(X1(end),Y1(end),Z1(end),'ro','MarkerFaceColor','r');
plot3(X2(end),Y2(end),Z2(end),'go','MarkerFaceColor','g');
axis equal;
grid on;
title('Motion of Bodies relative to the Inertial Reference Frame');
xlabel('Inertial X (km)');
ylabel('Inertial Y (km)');
zlabel('Inertial Z (km)');
%legend('Body1');
%% Motion of two bodies as seen from the Barycenter
figure(2)
hold on
plot3(X1-X_CM,Y1-Y_CM,Z1-Z_CM,'r');hold on % body 1
plot3(X2-X_CM,Y2-Y_CM,Z2-Z_CM,'g'); % body 2
grid on;
axis equal;
grid on;
title('Motion of Bodies relative to the Barycenter');
xlabel('Inertial X (km)');
ylabel('Inertial Y (km)');
zlabel('Inertial Z (km)');

view([X_CM(1) Y_CM(1) Z_CM(1)])
%view([0 0 0])
%% Animates the combined motion of both the masses
figure(3)
plot3(0,0,0,'ko','MarkerFaceColor','k');hold on % inertial origin
title('Motion of Bodies relative to the Inertial Reference Frame Animation');
xlabel('Inertial X (km)');
ylabel('Inertial Y (km)');
zlabel('Inertial Z (km)');
axis equal
for i = 2:Steps
%     calc = sqrt((X2(i)-X1(i))^2 + (Y2(i)-Y1(i))^2 + (Z2(i)-Z1(i))^2);
%     if calc < 200
%         break
%     end
    plot3(X1(i),Y1(i),Z1(i),'r');hold on % body 1
    plot3(X2(i),Y2(i),Z2(i),'g'); % body 2
    plot3([X1(i-1) X1(i)],[Y1(i-1) Y1(i)],[Z1(i-1) Z1(i)],'r');
    plot3([X2(i-1) X2(i)],[Y2(i-1) Y2(i)],[Z2(i-1) Z2(i)],'g');grid on
    pause(0.01)
    
    
end
%% RK 4 SCHEME
    function y = RK4(Derivative, t0, y0, h)
        % Runga Kutta to solve 2BEOM
        k1 = Derivative(t0    , y0          );
        k2 = Derivative(t0+h/2, y0+(h/2)*k1);
        k3 = Derivative(t0+h/2, y0+(h/2)*k2);
        k4 = Derivative(t0+h  , y0+    h*k3);
        y = y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
%% Differential equation in state space form ( Newton law of universal gravitation )
    function dydt = Derivative(~,y)
        r = sqrt((y(4)-y(1))^2 + (y(5)-y(2))^2 + (y(6)-y(3))^2); % distance between bodies 
        % dydt = [velocities and accelerations]
        dydt = [y(7),y(8),y(9),y(10),y(11),y(12),((G*m2)/(r^3))*(y(4)-y(1)),((G*m2)/(r^3))*(y(5)-y(2)),((G*m2)/(r^3))*(y(6)-y(3)),((G*m1)/(r^3))*(y(1)-y(4)),((G*m1)/(r^3))*(y(2)-y(5)),((G*m1)/(r^3))*(y(3)-y(6))]';
    end
end

