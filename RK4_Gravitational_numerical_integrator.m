function [X,Y,Z,VX,VY,VZ] = RK4_Gravitational_numerical_integrator(State_Vector,Time_span_hrs,Celestial_object_about_which_body_is_Orbiting)
format long;
if Celestial_object_about_which_body_is_Orbiting == "Sun"
    mu = 132712440041.93938; % km^3/s^2 sun gravitational parameter
elseif Celestial_object_about_which_body_is_Orbiting == "Earth"
    mu = 398600.435436; % km^3/s^2 earth gravitational parameter
end
y0 = State_Vector';% state vector of orbiting body
h=100;%step size for integration 
%initial values
t0 = 0; % sec
tend = Time_span_hrs*3600;% sec, iterating till the end of one time period
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
    y   = RK4(@deriv,t0, y0, h);
    y0 = y;
    t0 = t0+h;
    X(i) = y0(1);
    Y(i) = y0(2);
    Z(i) = y0(3);
    VX(i) = y0(4);
    VY(i) = y0(5);
    VZ(i) = y0(6);
end
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
    r=sqrt((y(1)^2)+(y(2)^2)+(y(3)^2));% calculated the L2 norm
    dydt=[y(4); y(5); y(6); -(mu/r^3)*y(1); -(mu/r^3)*y(2); -(mu/r^3)*y(3)];
end
end