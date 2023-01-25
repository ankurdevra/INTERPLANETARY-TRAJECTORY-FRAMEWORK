function [Right_Ascension,Declination] = RA_DEC_from_geocentric_equatiorial_position_vector(Position_vector)
% The following code calulates Right Ascension and Declination of an body
% orbiting earth using position vector measured in geocentric equatorial
% frame.
% REQUIRED INPUTS:
% Position_vector = Position vector of body as a [1X3] row vector [km]
% OUTPUTS:
% Right_Ascension = Right ascension of orbiting body in deg
% Declination = Declination of orbiting body in deg
%% Creator:- ANKUR DEVRA 
% Develope Date - 1 July 2022
% Iteration 1 -
%% INPUT:
X = Position_vector(1); % x component of position vector km
Y = Position_vector(2); % y component of position vector km
Z = Position_vector(3); % z component of position vector km
r = norm(Position_vector); % magnitude of poisition vector Km
%% CALCULATIONS:
% calculating direction cosines of position vector
l = X/r;
m = Y/r;
n = Z/r;
% calculating declination
delta = asind(n); % declination deg
% calculating right ascension
if m>0
    alpha = acosd(l/cosd(delta)); % right ascension deg
elseif m==0
    alpha = 360-acosd(l/cosd(delta)); % right ascension deg
elseif m<0
    alpha = 360-acosd(l/cosd(delta)); % right ascension deg
end
%% OUTPUT:
Right_Ascension = alpha; % deg
Declination = delta; %  deg
end