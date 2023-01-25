function [Universal_anomaly] = Universal_anomaly_Using_Newton_iteration(delta_t,initial_radius,initial_radial_velocity,semimajor_axis)
% The folling code calculates Universal anomaly using universal kepler
% equation utilizing stumpff functions.
% NEED MODIFICATION IN SEMI MAJOR AXIS FOR ALL OTHER ORBITS,
% CURRENTLY ONLY WORKS FOR HYPERBOLIC ORBITS.
% REQUIRED INPUTS:
% delta_t = Chnage in time sec.
% initial_radius = km, initial distance of secondary body 
% initial_radial_velocity = km/sec, initial speed of secondary body
%semimajor_axis = km, semimajor axis of given orbit ( positive values only ) 
% OUTPUT:
% Universal_anomaly = Universal anomaly after given delta_t, km^1/2
%% Creator:- ANKUR DEVRA 
% Develope Date - 30 June 2022
% Iteration 1 -
%% STARTING VALUE:
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% INPUTS:
r0 = initial_radius;%km initial distance of secondary body 
vr0 = initial_radial_velocity;%km/sec initial speed of secondary body
a=-semimajor_axis;%km, (negative for hyperbolic orbit, MODIFY FOR ALL OTHER ORBITS) negative semi-major axis for hyperbolic orbit taken in accordance with universal kepler equation
%% CALCULATIONS:
alpha = 1/a; % reciprocal of semi-major axis
Chi0 = sqrt(mu_earth)*abs(alpha)*delta_t;% initial guess of universal anomaly
matrix = []; % initializing empty matrix
while abs(((((r0*vr0*(Chi0^2)*Stumpff_C(alpha*(Chi0)^2))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^3)*Stumpff_S(alpha*(Chi0)^2) + r0*Chi0 - sqrt(mu_earth)*delta_t)/((((r0*vr0*Chi0*(1-alpha*((Chi0)^2)*Stumpff_S(alpha*(Chi0)^2))))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^2)*Stumpff_C(alpha*(Chi0)^2) + r0))) > 10^(-8) % specifying error tolerance
    ratio = ((((r0*vr0*(Chi0^2)*Stumpff_C(alpha*(Chi0)^2))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^3)*Stumpff_S(alpha*(Chi0)^2) + r0*Chi0 - sqrt(mu_earth)*delta_t)/((((r0*vr0*Chi0*(1-alpha*((Chi0)^2)*Stumpff_S(alpha*(Chi0)^2))))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^2)*Stumpff_C(alpha*(Chi0)^2) + r0));
    Chi0 = Chi0 - ratio; % Universal anomaly using newton iteration
    matrix = [matrix;Chi0];% stores the result in a matrix after each iteration
end
%% OUTPUT
Universal_anomaly = matrix(end,1);
% [row,~] = size(matrix);
% iterations = 1:row;
% Anomaly_table = [iterations' matrix];
% a2t = array2table(Anomaly_table,"VariableNames",["Iteration","Universal Anomaly (Km^(1/2))"]);disp(a2t)
%% STUMPFF FUNCTION 
    function S = Stumpff_S(z)
        if z>0
            S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
        elseif z<0
            S = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z))^3;
        elseif z==0
            S = 1/6;
        end
    end
    function C = Stumpff_C(z)
        if z>0
            C = (1-cos(sqrt(z)))/((z));
        elseif z<0
            C = (cosh(sqrt(-z))-1)/((-z));
        elseif z==0
            C = 1/2;
        end
    end

end