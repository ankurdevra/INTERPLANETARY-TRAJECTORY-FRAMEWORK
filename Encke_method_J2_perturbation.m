function [] = Encke_method_J2_perturbation(Initial_Perigee,Initial_Apogee,Initial_RAAN,Initial_Inclination,Initial_Argument_of_Perigee,Initial_True_Anomaly,Elapsed_time_hrs)
% The following code calculates the variation of orbitial parameter for an
% earth orbiting body due to J2 zonal perturbation.
% NEEDS A LOT OF WORK, VERY UNSTABLE AND UNRELIABLE CODE
% WORK REQUIRED
% MAY HAVE MULTIPLE CALCULATION ERRORS
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal constant
J2 = 0.00108263; % J2 zonal harmonics of earth accounting for gravitational perturbation 
R_earth_radius = 6378.137; % km Earth Equatorial Radius
%% CALCULATIONS
Initial_Eccentricity = (Initial_Apogee-Initial_Perigee)/(Initial_Apogee+Initial_Perigee); % current inclination of orbit around earth
Initial_semimajoraxis = (Initial_Apogee+Initial_Perigee)/2; % km, initial semimajor axis of body around earth
Initial_Angular_momentum = sqrt(mu_earth*Initial_semimajoraxis*(1-Initial_Eccentricity^2)); % km/sec^2, initial angular momentum of body around earth orbit
Initial_Time_Period = (2*pi)*sqrt(((Initial_semimajoraxis)^3)/mu_earth); % sec, initial time period of earth orbiting body
% Now calculating the initial Geocentric Equatorial state vectors of body around earth by
% invoking Geocentric_equi_and_perifocal_state_vector_from_orbital_element
% function
% these state vectors are at time t0 = 0
[~,Geocentric_equatorial_state_vectors] = Geocentric_equi_and_perifocal_state_vector_from_orbital_element(Initial_Angular_momentum,Initial_Eccentricity,Initial_Inclination,Initial_RAAN,Initial_Argument_of_Perigee,Initial_True_Anomaly);
r0 = Geocentric_equatorial_state_vectors(1:3);% km, [1X3] initial position vector of orbiting body around earth
v0 = Geocentric_equatorial_state_vectors(4:6);% km/sec, [1X3] initial velocity vector of orbiting body around earth

% Encke integration data
t0 = 0; % sec, initial starting time
tf = Elapsed_time_hrs*3600; % sec, final ending time
delta_T = Initial_Time_Period/10; % sec, encke method time step
y0 = [r0 v0]; % [1X6], km and km/sec, initial state vector to earth orbiting body
del_y0 = zeros(1,6); % [1X6], km and km/sec, initial perturbation state vector set to 0
t = t0;
t = t+delta_T; % sec, first time step for encke method
t_equi = [];state_vectors=[];
while t <= (tf+delta_T) % sec, time step
    [Time,Sol] = ode45(@perturbation,[t0 t],del_y0); % calls ODE45 to solve the perturbation differential equation from t0 to t with initial values del_y0
    % Now calculating osculating state vector at time t using initial state
    % vector y0
    [Osculating_state_vector,~] = future_state_vector_from_initial_using_universal_anomaly(y0,t-t0);
    R_osc = Osculating_state_vector(1:3); % km, [1X3] osculating position vector at time t
    V_osc = Osculating_state_vector(4:6); % km/sec,[1X3] osculating velocity vector at time t
    r0 = R_osc+Sol(end,1:3); % km, [1X3] rectification, updating the position vector of body
    v0 = V_osc+Sol(end,4:6); % km/sec, [1X3] rectification, updating the velocity vector of body
    y0 = [r0 v0]; % km and km/sec, [1X6] updated state vector of body after accounting for gravitational perturbation
    t0 = t; % sec, updating time
    t = t+delta_T; % sec, updated time step for next iteration
    t_equi = [t_equi;t]; % sec, matrix of equidistant time steps
    state_vectors = [state_vectors;y0]; % stores the state vectors at each time step
    
end

for i = 1:length(t_equi)
    [Angular_momentum(i),Inclination(i),Eccentricity(i),RAAN(i),Argument_of_Perigee(i),True_anomaly(i),~,~,~,~,Semimajor_axis(i),~,~,~,~] = Orbital_elements_from_State_vectors(state_vectors(i,1:6));
        
end
% PERTURBATION SUBFUCTION FOR ENCKE METHOD
    function dfdt = perturbation(~,f)
        del_r = f(1:3); % km, [1X3] deviation in position, initially 0
        del_v = f(4:6); % km/sec, [1X3] deviation in velocity, initially 0
        % Now calculating osculating state vector at time t using initial state
        % vector y0
        [Osculating_state_vector,~] = future_state_vector_from_initial_using_universal_anomaly(y0,t-t0);
        R_osc = Osculating_state_vector(1:3); % km, [1X3] osculating position vector at time t
        V_osc = Osculating_state_vector(4:6); % km/sec,[1X3] osculating velocity vector at time t
        R_per = R_osc+del_r; % km, perturbed (actual) position vector [1X3]
        r = norm(R_per); % km, magnitude of perturbed (actual) position vector
        V_per = V_osc+del_v; % km/sec, perturbed (actual) velocity vector [1X3]
        x = R_per(1);y = R_per(2);z = R_per(3); % km, x,y,z components of perturbed (actual) position vector
        p = (((3*J2*mu_earth*((R_earth_radius)^2))/(2*((r)^4)))).*[(x/r)*(((5*(z^2))/(r^2))-1) (y/r)*(((5*(z^2))/(r^2))-1) (z/r)*(((5*(z^2))/(r^2))-3)]; % [1X3] km/sec^2 perturbing gravitational acceleration vector due to J2
        del_a = (-mu_earth/((norm(R_osc))^3)).*(del_r - ((1-(((norm(R_osc))^3)/((r)^3))).*R_per)) + p; % km/sec^2, [1X3] total perturbing acceleration vector
        dfdt = [del_v(1) del_v(2) del_v(3) del_a(1) del_a(2) del_a(3)]'; % km/sec and km/sec^2 [6X1] derivative vector return to ODE
    end
end