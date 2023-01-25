clc;clear
[State_vectors_JPL_Horizons_Earth] = State_Vector_of_Planets_at_specified_Date_and_Time([2036 9 1],"Earth");
[State_vectors_JPL_Horizons_Jupiter] = State_Vector_of_Planets_at_specified_Date_and_Time([2036 9 1],"Jupiter");
[State_vectors_JPL_Horizons_Saturn] = State_Vector_of_Planets_at_specified_Date_and_Time([2036 9 1],"Saturn");
[State_vectors_JPL_Horizons_Uranus] = State_Vector_of_Planets_at_specified_Date_and_Time([2036 9 1],"Uranus");

State_Vector_Earth = State_vectors_JPL_Horizons_Earth;%[-5.265930693018012E+07,1.379224284333397E+08,-1.202417292883992E+04,-3.011257164725201E+01,-1.008747579302776E+01,2.500762775889611E-01];% state vector of earth
State_Vector_Jupiter = State_vectors_JPL_Horizons_Jupiter;%[5.932409994340167E+08,4.482020205876219E+08,-1.518150458725777E+07,-9.307071935861673E+00,1.258725083995466E+01,3.371409490836612E-01]; % state vector of jupiter
State_Vector_Saturn = State_vectors_JPL_Horizons_Saturn;
State_Vector_Uranus = State_vectors_JPL_Horizons_Uranus;
% ORBITAL ELEMENT DATA
% [Angular_momentum_Earth,Inclination_Earth,Eccentricity_Earth,RAAN_Earth,Argument_of_Perigee_Earth,True_anomaly_Earth,Radial_velocity_Earth,Time_period_of_orbit_hrs_Earth,Radius_perigee_Earth,Radius_apogee_Earth,Semimajor_axis_Earth,Semilatus_rectum_Earth,Eccentric_anomaly_Earth,Mean_anomaly_Earth,Orbits_in_a_day_Earth] = Heliocentric_Orbital_elements_from_State_vectors(State_Vector_Earth)
% [Angular_momentum_Jupiter,Inclination_Jupiter,Eccentricity_Jupiter,RAAN_Jupiter,Argument_of_Perigee_Jupiter,True_anomaly_Jupiter,Radial_velocity_Jupiter,Time_period_of_orbit_hrs_Jupiter,Radius_perigee_Jupiter,Radius_apogee_Jupiter,Semimajor_axis_Jupiter,Semilatus_rectum_Jupiter,Eccentric_anomaly_Jupiter,Mean_anomaly_Jupiter,Orbits_in_a_day_Jupiter] = Heliocentric_Orbital_elements_from_State_vectors(State_Vector_Jupiter)
% [Angular_momentum_Saturn,Inclination_Saturn,Eccentricity_Saturn,RAAN_Saturn,Argument_of_Perigee_Saturn,True_anomaly_Saturn,Radial_velocity_Saturn,Time_period_of_orbit_hrs_Saturn,Radius_perigee_Saturn,Radius_apogee_Saturn,Semimajor_axis_Saturn,Semilatus_rectum_Saturn,Eccentric_anomaly_Saturn,Mean_anomaly_Saturn,Orbits_in_a_day_Saturn] = Heliocentric_Orbital_elements_from_State_vectors(State_Vector_Saturn)
% [Angular_momentum_Uranus,Inclination_Uranus,Eccentricity_Uranus,RAAN_Uranus,Argument_of_Perigee_Uranus,True_anomaly_Uranus,Radial_velocity_Uranus,Time_period_of_orbit_hrs_Uranus,Radius_perigee_Uranus,Radius_apogee_Uranus,Semimajor_axis_Uranus,Semilatus_rectum_Uranus,Eccentric_anomaly_Uranus,Mean_anomaly_Uranus,Orbits_in_a_day_Uranus] = Heliocentric_Orbital_elements_from_State_vectors(State_Vector_Uranus)
%
Time_span_hrs_Earth=365*24;
% Time_span_hrs_Jupiter = Time_span_hrs_Earth*20;
% Time_span_hrs_Saturn = Time_span_hrs_Earth*40;
% Time_span_hrs_Uranus = Time_span_hrs_Earth*60;

Celestial_object_about_which_body_is_Orbiting = "Sun";
[X_Earth,Y_Earth,Z_Earth,VX_Earth,VY_Earth,VZ_Earth] = RK4_Gravitational_numerical_integrator(State_Vector_Earth,Time_span_hrs_Earth,Celestial_object_about_which_body_is_Orbiting);
% [X_Jupiter,Y_Jupiter,Z_Jupiter,VX_Jupiter,VY_Jupiter,VZ_Jupiter] = RK4_Gravitational_numerical_integrator(State_Vector_Jupiter,Time_span_hrs_Jupiter,Celestial_object_about_which_body_is_Orbiting);
% [X_Saturn,Y_Saturn,Z_Saturn,VX_Saturn,VY_Saturn,VZ_Saturn] = RK4_Gravitational_numerical_integrator(State_Vector_Saturn,Time_span_hrs_Saturn,Celestial_object_about_which_body_is_Orbiting);
% [X_Uranus,Y_Uranus,Z_Uranus,VX_Uranus,VY_Uranus,VZ_Uranus] = RK4_Gravitational_numerical_integrator(State_Vector_Uranus,Time_span_hrs_Uranus,Celestial_object_about_which_body_is_Orbiting);

%% plotting SUN
imData = imread('sun_mercator_projection.jpg');
[xS,yS,zS] = sphere(50);
sun_radius = (6.6846e-9)*696340;  % AU
xSE = sun_radius*xS;
ySE = sun_radius*yS;
zSE = sun_radius*zS;
figure(1)
surface(xSE,ySE,zSE);hold on
axis equal;
ch = get(gca,'children');
set(ch,'facecolor','texturemap','cdata',flipud(imData),'edgecolor','black')
set(gca,'color','white')

%% plotting planetary orbit in AU
figure(2)
plot3((6.6846e-9).*X_Earth,(6.6846e-9).*Y_Earth,(6.6846e-9).*Z_Earth,'b','LineWidth',2);hold on
plot3((6.6846e-9).*X_Jupiter,(6.6846e-9).*Y_Jupiter,(6.6846e-9).*Z_Jupiter,'r','LineWidth',2);
plot3((6.6846e-9).*X_Saturn,(6.6846e-9).*Y_Saturn,(6.6846e-9).*Z_Saturn,'k','LineWidth',2);
plot3((6.6846e-9).*X_Uranus,(6.6846e-9).*Y_Uranus,(6.6846e-9).*Z_Uranus,'b','LineWidth',2);

title('Planetary Orbits Reference frame : Ecliptic of J2000.0');
xlabel('X (AU)')
ylabel('Y (AU)')
zlabel('Z (AU)')
grid on ;
axis vis3d
legend('Earth','Jupiter','Saturn','Uranus')
