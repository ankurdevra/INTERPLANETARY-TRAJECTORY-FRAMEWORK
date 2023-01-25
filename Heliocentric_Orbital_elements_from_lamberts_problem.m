function [Angular_momentum,Inclination,Eccentricity,RAAN,Argument_of_Perigee,True_anomaly,Radial_velocity,Time_period_of_orbit_hrs,Radius_perigee,Radius_apogee,Semimajor_axis,Semilatus_rectum,Eccentric_anomaly,Mean_anomaly,Orbits_in_a_day,state_vector_at_point_1,state_vector_at_point_2] = Heliocentric_Orbital_elements_from_lamberts_problem(r1_vec,r2_vec,elapsed_time_hrs,Orbit_type)
% ANKUR DEVRA MAE 460
% The following code calculates the orbital elements of a body around sun
% from two successive heliocentric position vectors who have a elapsed time
% in between them by solving the lamberts problem.
% REQUIRED INPUTS:
% r1_vec = km, heliocentric position vector for first observation
% r2_vec = km, heliocentric position vector for second observation
% elapsed_time_hrs = hrs, elapsed time between two geocentric position
% observations
% Orbit_type = 'prograde' OR 'retrograde' orbit type around earth
% OUTPUTS:
% Angular_momentum = specific angular momentum km^2/sec
% Inclination = inclination of orbit is degrees, if i>90 retrograde orbit
% Eccentricity =  eccentricity
% RAAN =  Right Ascension of Ascending Node, degrees
% Argument_of_Perigee =  argument of perigee degrees
% True_anomaly =  true anomaly degrees
% Radial_velocity =  km/sec (v_r) if v_r>0 trajectory away from preigee if v_r<0 towards perigee
% Time_period_of_orbit_hrs = Time period of orbit in hrs
% Radius_perigee =  radius of perigee in km
% Radius_apogee =  radius of apogee in km
% Semimajor_axis =  semi-major axis in km
% Semilatus_rectum =  semi-latus rectum in km
% Eccentric_anomaly =  Eccentric anomaly in degrees
% Mean_anomaly =  Mean anomaly in degrees
% Orbits_in_a_day =  number of orbits in a day 
%% Creator:- ANKUR DEVRA 
% Develope Date - 3 July 2022
% Iteration 1 -
%% Starting data
mu_sun = 132712440041.93938;%1.32712440018*10^(11); % km^3/s^2 sun gravitational parameter
%% CALCULATIONS
elapsed_time_sec = elapsed_time_hrs*3600; % sec, time elapsed between two successive position meausrement in geocentric position
r1_mag = norm(r1_vec); % km magnitude of r1 position vector
r2_mag = norm(r2_vec); % km magnitude of r2 position vector

position_vectors_cross_product = cross(r1_vec,r2_vec); % km^2 cross product of given position vectors
z_component = position_vectors_cross_product(3); % km^2, z component of position vector cross product

if Orbit_type == "prograde"
    if z_component>=0
        delta_theta = acosd((dot(r1_vec,r2_vec))/(r1_mag*r2_mag)); % deg change in true anomaly
    elseif z_component<0
        delta_theta = 360-acosd((dot(r1_vec,r2_vec))/(r1_mag*r2_mag)); % deg change in true anomaly
    end
elseif Orbit_type == "retrograde"
    if z_component<0
        delta_theta = acosd((dot(r1_vec,r2_vec))/(r1_mag*r2_mag)); % deg change in true anomaly
    elseif z_component>=0
        delta_theta = 360-acosd((dot(r1_vec,r2_vec))/(r1_mag*r2_mag)); % deg change in true anomaly
    end
end

A = sind(delta_theta)*sqrt((r1_mag*r2_mag)/(1-cosd(delta_theta))); % km intermediate constant in solving lamberts problem
z = -100; % starting guess value
% checking with while loop when F(z) approx changes sign, this will be the
% starting value of z for newton iteration method
while F(z)<0
    z = z+0.1;
end
% Doing newton iteration to solve for z. sign of z determines the orbit
% type. z<0 hyperbola; z=0 parabola; z>0 ellipse
while abs(F(z)/dFdz(z)) > 10^(-8) % specifying error tolerance
    z = z - (F(z)/dFdz(z));
end

% Calculating lagrange coefficients
f = 1-(y(z)/r1_mag); % f lagrange coeff
g = A*sqrt(y(z)/mu_sun); % sec, g lagrange coeff
f_dot = ((sqrt((mu_sun))/(r1_mag*r2_mag))*(sqrt(y(z)/Stumpff_C(z))))*((z*Stumpff_S(z))-1); % sec, f_dot lagrange coeff
g_dot = 1-(y(z)/r2_mag); % g_dot lagrange coeff

% Calculating state vectors at point 1 and point 2 using lagrange
% coefficients

v1_vec = (1/g).*(r2_vec - (f.*(r1_vec))); % km/sec velocity vector at point 1
v2_vec = (1/g).*((g_dot.*(r2_vec)) - r1_vec); % km/sec velocity vector at point 2

state_vector_at_point_1 = [r1_vec,v1_vec]; % KM AND KM/SEC heliocentric state vector at point 1 (start of observation)
state_vector_at_point_2 = [r2_vec,v2_vec]; % KM AND KM/SEC heliocentric state vector at point 2 (end of observation)

%% OUTPUT
% invokes the Orbital_elements_from_State_vectors function to find the
% orbital elements of given orbit around sun.
% the anomalies are with respect to first observation
[Angular_momentum,Inclination,Eccentricity,RAAN,Argument_of_Perigee,True_anomaly,Radial_velocity,Time_period_of_orbit_hrs,Radius_perigee,Radius_apogee,Semimajor_axis,Semilatus_rectum,Eccentric_anomaly,Mean_anomaly,Orbits_in_a_day] = Heliocentric_Orbital_elements_from_State_vectors(state_vector_at_point_1);

%% SUBFUNCTIONS for lamberts problem
    function y_z = y(z)
        y_z = r1_mag + r2_mag + A*(((z*Stumpff_S(z))-1)/(sqrt(Stumpff_C(z)))); % km intermediate function in solving lamberts problem
    end
    function F_z = F(z)
        F_z = (((y(z))/(Stumpff_C(z)))^(3/2))*Stumpff_S(z) + A*sqrt(y(z)) - (sqrt(mu_sun))*elapsed_time_sec; % F(z) function of lamberts problem
    end

    function F_dash_z = dFdz(z)
        if z==0
            F_dash_z = ((sqrt(2)/40)*(y(0))^(3/2)) + (A/8)*((sqrt(y(0))) + A*(1/(2*y(0))));
        else
            F_dash_z = (((y(z))/(Stumpff_C(z)))^(3/2))*(((1/(2*z))*(Stumpff_C(z)-((3*(Stumpff_S(z))))/(2*Stumpff_C(z)))) + (((3*(Stumpff_S(z))^2)/(4*Stumpff_C(z)))))...
                + (A/8)*((((3*Stumpff_S(z))/(Stumpff_C(z)))*sqrt(y(z))) + A*sqrt((Stumpff_C(z))/(y(z))));
        end
    end
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