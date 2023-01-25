function [delta_v_at_A,delta_v_at_B,total_delta_v,Total_propellant_required,TOF_hohmann] = Hohmann_tranfer_from_lower_to_higher_orbit(mass_of_spacecraft,initial_orbit_perigee,initial_orbit_apogee,final_orbit_perigee,final_orbit_apogee,specific_impulse_of_propellant)
% The following code calculates the the total delta-v and propellant
% required to execute a hohmann transfer from a lower orbit 1 to final
% circular orbit 3.
% apogee and perigee distance are wrt earth's center
% REQUIRED INPUTS:
% mass_of_spacecraft = kg, initial mass of spacecraft (body+fuel)
% initial_orbit_perigee = km, perigee of initial orbit of body around earth
% initial_orbit_apogee = km, apogee of initial orbit of body around earth
% final_orbit_perigee = km, perigee of final orbit of body around earth
% final_orbit_apogee = km, apogee of final orbit of body around earth
% specific_impulse_of_propellant = sec, specific impulse of propellant
% final orbit must be a cicular orbit
% OUTPUTS:
% delta_v_at_A = km/sec delta-v required at perigee of orbit 1 to get into an elliptical trajectory of transfer orbit 2
% delta_v_at_B = km/sec (apogee kick) delta-v required at point B (radius of final circular orbit 3) to get into the final circular orbit 3 from tranfer elliptical orbit 2
% total_delta_v = km/sec, total delta-v required to go from lower orbit 1 to final circular orbit 3
% Total_propellant_required = kg change in mass of spacecraft after buring fuel (total propellant required) to execute hohmann tranfer
% TOF_hohmann = sec, time of flight in hohmann trajectory
%% Creator:- ANKUR DEVRA 
% Develope Date - 30 June 2022
% Iteration 1 -
%% Starting data
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal parameter
g0 = (9.81)/1000; % km/sec^2 sea level standard graviational accelaration 
%% INPUTS
m = mass_of_spacecraft; % kg initial mass of spacecraft (body+fuel)
rp1 = initial_orbit_perigee; % km, perigee of initial orbit of body around earth
ra1 = initial_orbit_apogee; % KM, apogee of initial orbit of body around earth
rp2 = final_orbit_perigee; % km, perigee of final orbit of body around earth
ra2 = final_orbit_apogee; % km, apogee of final orbit of body around earth
Isp = specific_impulse_of_propellant; % sec, specific impulse of propellant
%% CALCULATIONS
% initial orbit 1
h1 = sqrt(2*mu_earth)*(sqrt((ra1*rp1)/(ra1+rp1))); % km/sec^2; angular momentum of initial orbit 1
% transfer orbit 2
% The tranfer orbit will have same perigee as initial orbit 1 but different
% apogee. Tranfer orbit apogee will be the radius of final circular orbit
h2 = sqrt(2*mu_earth)*(sqrt((ra2*rp1)/(ra2+rp1))); % km/sec^2; angular momentum of transfer elliptical orbit 
% final orbit 3
h3 = sqrt(mu_earth*rp2); % km/sec^2, angular momentum of final circular orbit 3

va_1 = h1/rp1; % km/sec speed on orbit 1 at point A (perigee of orbit 1)
va_2 = h2/rp1; % km/sec speed on tranfer orbit 2 at point A (perigee of elliptical tranfer orbit 2)
delta_v_a = va_2-va_1; % km/sec delta-v required at perigee of orbit 1 to get into an elliptical trajectory of transfer orbit 2

v_b_2 = h2/rp2; % km/sec speed on elliptical orbit 2 at point B (radius of final circular orbit 3)
v_b_3 = h3/rp2; % km/sec speed on circular orbit 3 at point B (radius of final circular orbit 3)
delta_v_b = v_b_3-v_b_2; % km/sec (apogee kick) delta-v required at point B (radius of final circular orbit 3) to get into the final circular orbit 3 from tranfer elliptical orbit 2

total_delta_v = abs(delta_v_a)+abs(delta_v_b); % km/sec, total delta-v required to go from lower orbit 1 to final circular orbit 3

delta_m = m*(1 - exp(-total_delta_v/(Isp*g0))); % kg change in mass of spacecraft after buring fuel (total propellant required) to execute hohmann tranfer  
a = (1/2)*(rp1 + ra2); % km, semimajor axis of hohmann ellipse
%% OUTPUT
delta_v_at_A = delta_v_a; % km/sec delta-v required at perigee of orbit 1 to get into an elliptical trajectory of transfer orbit 2
delta_v_at_B = delta_v_b; % km/sec (apogee kick) delta-v required at point B (radius of final circular orbit 3) to get into the final circular orbit 3 from tranfer elliptical orbit 2
total_delta_v = total_delta_v; % km/sec, total delta-v required to go from lower orbit 1 to final circular orbit 3
Total_propellant_required = delta_m; % kg change in mass of spacecraft after buring fuel (total propellant required) to execute hohmann tranfer 
TOF_hohmann = (1/2)*(((2*pi)/sqrt(mu_earth))*a^(3/2)); % sec, time of flight in hohmann trajectory
end