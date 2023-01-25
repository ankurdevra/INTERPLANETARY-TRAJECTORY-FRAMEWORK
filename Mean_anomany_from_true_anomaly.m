function [Mean_anomaly] = Mean_anomany_from_true_anomaly(True_anomaly,Orbital_eccentricity)
% The following code converts given true anomaly of an orbit to its
% corresponing mean anomaly.
% REQUIRED INPUTS:
% True_anomaly = True anomaly which we wish to convert to mean anomaly in radians.
% Orbital_eccentricity = eccentricty of given orbit.
% OUTPUT:
% Mean_anomaly = Mean anomaly at given true anomaly in radians 
%% Creator:- ANKUR DEVRA 
% Develope Date - 23 May 2022
% Iteration 1 -
%% INPUTS:
theta = True_anomaly;% radian, true anomaly 
e = Orbital_eccentricity;% eccentricity of orbit
%% OUTPUT
Mean_anomaly = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2)) - (e*sqrt(1-(e)^2)*sin(theta))/(1+e*cos(theta));
Me = ['The mean anomaly is ',num2str(Mean_anomaly),' radians'];disp(Me);

end
