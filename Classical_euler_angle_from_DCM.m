function [alpha,beta,gamma] = Classical_euler_angle_from_DCM(DCM)
% The following code is used to find the classical euler angles form
% symmetric classical euler angle DCM
% * MAY HAVE tan QUADRANT ANBIGUITY NEEDS CORRECTION
% REQUIRED INPUTS:
% DCM = [Q11 Q12 Q13 Q21 Q22 Q23 Q31 Q32 Q33] enter as a row matrix, row
% wise
% OUTPUTS:
% alpha = first rotation angle about 3rd axis deg
% beta = second rotation angle about 1st axis deg
% gamma = third rotation angle about 3rd axis deg
%% Creator:- ANKUR DEVRA 
% Develope Date - 2 July 2022
% Iteration 1 -
%% CALCULATIONS
% classical euler sequence DCM
DCM = [DCM(1,1) DCM(1,2) DCM(1,3);DCM(1,4) DCM(1,5) DCM(1,6);DCM(1,7) DCM(1,8) DCM(1,9)];
Q31 = DCM(3,1);
Q32 = DCM(3,2);
Q33 = DCM(3,3);
Q13 = DCM(1,3);
Q23 = DCM(2,3);
%% OUTPUT
alpha = atand(Q31/-Q32); % first rotation angle about 3rd axis deg
if alpha<0
    alpha = alpha + 360; % first rotation angle about 3rd axis deg
end
beta = acosd(Q33); % second rotation angle about 1st axis deg
gamma = atand(Q13/Q23); % third rotation angle about 3rd axis deg
if gamma<0
    gamma = gamma + 360; % third rotation angle about 3rd axis deg
end
end