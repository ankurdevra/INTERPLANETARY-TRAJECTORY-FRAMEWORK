function [future_state_vector,future_position_and_velocity] = future_state_vector_from_initial_using_universal_anomaly(initial_state_vector,delta_t)
% The following code calculates future state vector of a body orbiting
% earth from given initial state vectors after elapsed time. Also
% calculates the final position and velocity of body after elapsed time.
% It uses Universal keplers equation and universal anomaly by utilising
% lagrange coefficients using Newton method to iteratively solve for universal anomaly using stumpff
% functions
% WORKS FOR ANY TYPE OF ORBIT.
% REQUIRED INPUTS:
% initial_state_vector = initial state vector to be input as [1X6] row matrix km and km/sec
% delta_t = elapsed time or change in time sec.
% OUTPUT:
% future_state_vector =  future state vector after elapsed time as [1X6] row matrix km and km/sec
% future_position_and_velocity = future position and velocity after elapsed time as [1X2] row vector km anf km/sec
%% Creator:- ANKUR DEVRA 
% Develope Date - 30 June 2022
% Iteration 1 -
%% STARTING VALUE:
format long
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
%% INPUTS:
initial_position_vector = initial_state_vector(1:3); % intial position vector km
initial_velocity_vector = initial_state_vector(4:6); % intial velocity vector km/sec
%% CALCULATIONS:
initial_position = norm(initial_position_vector); % magnitude of initial position vector km
initial_velocity = norm(initial_velocity_vector); % magnitude of initial velocity vector km/sec
initial_radial_velocity = (dot(initial_velocity_vector,initial_position_vector))/(initial_position); % initial radial velocity Km/sec
vr0 = initial_radial_velocity; % initial radial velocity Km/sec
r0 = initial_position; % magnitude of initial position vector km
alpha = (2/initial_position) - ((initial_velocity)^2)/mu_earth; % km^-1 determine orbit type. alpha>0 ellipse; alpha = 0 parabolic; alpha<0 hyperbola

% Newton method to iteratively solve for universal anomaly using stumpff
% function and universal kepler equation.
Chi0 = sqrt(mu_earth)*abs(alpha)*delta_t;% initial guess of universal anomaly
matrix = []; % initializing empty matrix
while abs(((((r0*vr0*(Chi0^2)*Stumpff_C(alpha*(Chi0)^2))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^3)*Stumpff_S(alpha*(Chi0)^2) + r0*Chi0 - sqrt(mu_earth)*delta_t)/((((r0*vr0*Chi0*(1-alpha*((Chi0)^2)*Stumpff_S(alpha*(Chi0)^2))))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^2)*Stumpff_C(alpha*(Chi0)^2) + r0))) > 10^(-8) % specifying error tolerance
    ratio = ((((r0*vr0*(Chi0^2)*Stumpff_C(alpha*(Chi0)^2))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^3)*Stumpff_S(alpha*(Chi0)^2) + r0*Chi0 - sqrt(mu_earth)*delta_t)/((((r0*vr0*Chi0*(1-alpha*((Chi0)^2)*Stumpff_S(alpha*(Chi0)^2))))/(sqrt(mu_earth))) + (1-(alpha*r0))*((Chi0)^2)*Stumpff_C(alpha*(Chi0)^2) + r0));
    Chi0 = Chi0 - ratio; % Universal anomaly using newton iteration
    matrix = [matrix;Chi0];% stores the result in a matrix after each iteration
end
Universal_anomaly = matrix(end,1); % Universal anomaly km^1/2
Chi = Universal_anomaly; % Universal anomaly km^1/2
z = alpha*Chi^2; % final argument of stumpff function to put in lagrange coeffcients.
% lagrange coefficients
f = 1 - (((Chi)^2)/r0)*Stumpff_C(z);
g = delta_t - (1/sqrt(mu_earth))*((Chi)^3)*Stumpff_S(z); % sec
final_position_vector = f.*initial_position_vector + g.*initial_velocity_vector; % km, final position vector
final_position = norm(final_position_vector); % magnitude of final position vector km
r = final_position; % final position km
f_dot = ((sqrt(mu_earth))/(r*r0))*(alpha*((Chi)^3)*Stumpff_S(z)-Chi); % sec^-1
g_dot = 1 - (((Chi)^2)/r)*Stumpff_C(z);
final_velocity_vector = f_dot.*initial_position_vector + g_dot.*initial_velocity_vector; % km.sec, final velocity vector
final_velocity = norm(final_velocity_vector); % magnitude of final velocity vector km/sec
v = final_velocity; % final velocity km/sec 
%% OUTPUT:
future_state_vector = [final_position_vector,final_velocity_vector]; % final state row (1X6) vector matrix km and km/sec
future_position_and_velocity = [r,v]; % final position and vector after elapsed time km and km/sec
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
