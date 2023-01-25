%% ANKUR DEVRA
% Solves for mean anomaly and true anomaly using newton iteration method to
% solve keplers equation
% and displays the result in a table.
% MIGHT NEED TO RECTIFY THE RANGE OF TRUE ANOMALY. LINE 19
function [eccentric_anomaly,true_anomaly] = Eccentric_and_true_anomaly_from_mean_anomaly(Mean_anomaly,Eccentricity)
M = Mean_anomaly;% RADIAN, mean anomaly
e = Eccentricity; % eccentricity
if M > pi
    E = M - e/2; % initial guess
elseif M < pi
    E = M + e/2; % initial guess
end
%E = M;
matrix=[]; % initailize empty matrix
while abs(((E-e*sin(E)-M)/(1-e*cos(E)))) > 10^(-8) % specifying error tolerance
    E = E - ((E-e*sin(E)-M)/(1-e*cos(E))); % ECCENTRIC ANOMALY USING NEWTON iteration METHOD
    theta = 2*atan((sqrt((1+e)/(1-e)))*tan(E/2));%+360; % true anomaly for each value of eccentric anomaly
    matrix = [matrix;[E,theta+2*pi]];% stores the result in a matrix after each iteration
end
eccentric_anomaly = matrix(end,1);%rad
true_anomaly = matrix(end,2);%rad
% [row,column] = size(matrix);
% iterations = 1:row;
% Anomaly_table = [iterations' matrix];
% a2t = array2table(Anomaly_table,"VariableNames",["Iteration","Eccentric Anomaly","True Anomaly"]);disp(a2t)
end
