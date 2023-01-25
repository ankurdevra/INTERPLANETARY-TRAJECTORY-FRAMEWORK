function [Hyperbolic_eccentric_anomaly] = Hyperbolic_eccentric_anomaly_from_Hyperbolic_mean_anomaly(eccentricity,Hyperbolic_mean_anomaly)
% The following code converts given hyperbolic mean anomaly and eccentricity 
% of an hyperbolic orbit to its corresponing hyperbolic eccentric anomaly
% using newton iterative method by solving kepler equation for hyperbola.
% Can also convert hyperbolic eccentric anomaly to true anomaly
% REQUIRED INPUTS:
% Hyperbolic_mean_anomaly = Hyperbolic mean anomaly which we wish to convert to hyperbolic eccentric anomaly in radians.
% eccentricity = eccentricty of given hyperbolic orbit.
% OUTPUT:
% Hyperbolic_eccentric_anomaly = hyperbolic eccentric anomaly at given hyperbolic mean anomaly and orbit eccentricity in radians 
%% Creator:- ANKUR DEVRA 
% Develope Date - 30 June 2022
% Iteration 1 -
%% INPUTS:
Mh = Hyperbolic_mean_anomaly;% radian, hyperbolic mean anomaly 
e = eccentricity;% eccentricity of  hyperbolic orbit
%% CALCULATIONS:
F0 = Mh; % initial guess of Hyperbolic_eccentric_anomaly
matrix=[]; % initailize empty matrix
while abs(((e*sinh(F0)-F0-Mh)/(e*cosh(F0)-1))) > 10^(-8) % specifying error tolerance
    F0 = F0 - ((e*sinh(F0)-F0-Mh)/(e*cosh(F0)-1)); % hyperbolic ECCENTRIC ANOMALY USING NEWTON iteration METHOD
    theta = 2*atan((sqrt((e+1)/(e-1)))*tanh(F0/2));%+360; % true anomaly for each value of hyperbolic eccentric anomaly
    matrix = [matrix;[F0,theta]];% stores the result in a matrix after each iteration
end
%% OUTPUT:
Hyperbolic_eccentric_anomaly = matrix(end,1);
% Hyperbolic_eccentric_anomaly = matrix(:,1);
% true_anomaly = matrix(:,2);
% [row,column] = size(matrix);
% iterations = 1:row;
% Anomaly_table = [iterations' matrix];
% a2t = array2table(Anomaly_table,"VariableNames",["Iteration","Hyperbolic Eccentric Anomaly (rad)","True Anomaly (rad)"]);disp(a2t)
end