clc;clear
Years = 2038:1:2040;
Months = 1:12;
diff = (Years(end)-Years(1))+1;
X_axis_vals = linspace(Years(1),Years(end)+1,diff*length(Months));
Date_Jan_1 = [2036 1 1 00 00 00];A = Date_Jan_1;
Date_Jan_15 = [2036 1 15 00 00 00];B = Date_Jan_15;
Date_Feb_1 = [2036 2 1 00 00 00];C = Date_Feb_1;
Date_Feb_15 = [2036 2 15 00 00 00];D = Date_Feb_15;
Date_Mar_1 = [2036 3 1 00 00 00];E = Date_Mar_1;
Date_Mar_15 = [2036 3 15 00 00 00];F = Date_Mar_15;
Date_Apr_1 = [2036 4 1 00 00 00];G = Date_Apr_1;
Date_Apr_15 = [2036 4 15 00 00 00];H = Date_Apr_15;
Date_May_1 = [2036 5 1 00 00 00];I = Date_May_1;
Date_May_15 = [2036 5 15 00 00 00];J = Date_May_15;
Date_June_1 = [2036 6 1 00 00 00];K = Date_June_1;
Date_June_15 = [2036 6 15 00 00 00];L = Date_June_15;
Date_July_1 = [2036 7 1 00 00 00];M = Date_July_1;
Date_July_15 = [2036 7 15 00 00 00];N = Date_July_15;
Date_Aug_1 = [2036 8 1 00 00 00];O = Date_Aug_1;
Date_Aug_15 = [2036 8 15 00 00 00];P = Date_Aug_15;
Date_Sept_1 = [2036 9 1 00 00 00];Q = Date_Sept_1;
Date_Sept_15 = [2036 9 15 00 00 00];R = Date_Sept_15;
Date_Oct_1 = [2036 10 1 00 00 00];S = Date_Oct_1;
Date_Oct_15 = [2036 10 15 00 00 00];T = Date_Oct_15;
Date_Nov_1 = [2036 11 1 00 00 00];U = Date_Nov_1;
Date_Nov_15 = [2036 11 15 00 00 00];V = Date_Nov_15;
Date_Dec_1 = [2036 12 1 00 00 00];W = Date_Dec_1;
Date_Dec_15 = [2036 12 15 00 00 00];X = Date_Dec_15;
Dates = [A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X];
%%
% for z = 1:length(Dates)/6
%     b=6*z;
%     kl = Dates((b-5):b)    
% end
%%
departure_hyperbolic_excess_velocity_Date_Jan_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Jan_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Jan_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Jan_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Feb_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Feb_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Feb_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Feb_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Mar_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Mar_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Mar_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Mar_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Apr_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Apr_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Apr_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Apr_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_May_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_May_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_May_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_May_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_June_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_June_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_June_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_June_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_July_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_July_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_July_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_July_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Aug_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Aug_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Aug_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Aug_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Sept_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Sept_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Sept_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Sept_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Oct_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Oct_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Oct_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Oct_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Nov_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Nov_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Nov_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Nov_15 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Dec_1 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Dec_1 = zeros(1,length(M));

departure_hyperbolic_excess_velocity_Date_Dec_15 = zeros(1,length(M));
arrival_hyperbolic_excess_velocity_Date_Dec_15 = zeros(1,length(M));
%% case 1: launch 1 Jan 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Jan_1(i,j),arrival_hyperbolic_excess_velocity_Date_Jan_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[A],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Jan_1 = departure_hyperbolic_excess_velocity_Date_Jan_1;D_Date_Jan_1 = Departure_Date_Jan_1';D_Date_Jan_1 = (D_Date_Jan_1(:))';
Arrival_Date_Jan_1 = arrival_hyperbolic_excess_velocity_Date_Jan_1;A_Date_Jan_1 = Arrival_Date_Jan_1';A_Date_Jan_1 = (A_Date_Jan_1(:))';
%figure(1)
%plot(X_axis_vals,D_Date_Jan_1);hold on;grid on;plot(X_axis_vals,A_Date_Jan_1); % case 1
%% case 2: launch 15 Jan 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Jan_15(i,j),arrival_hyperbolic_excess_velocity_Date_Jan_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[B],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Jan_15 = departure_hyperbolic_excess_velocity_Date_Jan_15;D_Date_Jan_15 = Departure_Date_Jan_15';D_Date_Jan_15 = (D_Date_Jan_15(:))';
Arrival_Date_Jan_15 = arrival_hyperbolic_excess_velocity_Date_Jan_15;A_Date_Jan_15 = Arrival_Date_Jan_15';A_Date_Jan_15 = (A_Date_Jan_15(:))';
%plot(X_axis_vals,D_Date_Jan_15);plot(X_axis_vals,A_Date_Jan_15); % case 2
%% case 3: launch 1 Feb 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Feb_1(i,j),arrival_hyperbolic_excess_velocity_Date_Feb_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[C],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Feb_1 = departure_hyperbolic_excess_velocity_Date_Feb_1;D_Date_Feb_1 = Departure_Date_Feb_1';D_Date_Feb_1 = (D_Date_Feb_1(:))';
Arrival_Date_Feb_1 = arrival_hyperbolic_excess_velocity_Date_Feb_1;A_Date_Feb_1 = Arrival_Date_Feb_1';A_Date_Feb_1 = (A_Date_Feb_1(:))';
%plot(X_axis_vals,D_Date_Feb_1);plot(X_axis_vals,A_Date_Feb_1); % case 3
%% Case 4: launch 15 feb 3026
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Feb_15(i,j),arrival_hyperbolic_excess_velocity_Date_Feb_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[D],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Feb_15 = departure_hyperbolic_excess_velocity_Date_Feb_15;D_Date_Feb_15 = Departure_Date_Feb_15';D_Date_Feb_15 = (D_Date_Feb_15(:))';
Arrival_Date_Feb_15 = arrival_hyperbolic_excess_velocity_Date_Feb_15;A_Date_Feb_15 = Arrival_Date_Feb_15';A_Date_Feb_15 = (A_Date_Feb_15(:))';
% plot(X_axis_vals,D_Date_Feb_15);plot(X_axis_vals,A_Date_Feb_15); % case 4
%% case 5: launch 1 march 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Mar_1(i,j),arrival_hyperbolic_excess_velocity_Date_Mar_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[E],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Mar_1 = departure_hyperbolic_excess_velocity_Date_Mar_1;D_Date_Mar_1 = Departure_Date_Mar_1';D_Date_Mar_1 = (D_Date_Mar_1(:))';
Arrival_Date_Mar_1 = arrival_hyperbolic_excess_velocity_Date_Mar_1;A_Date_Mar_1 = Arrival_Date_Mar_1';A_Date_Mar_1 = (A_Date_Mar_1(:))';
% plot(X_axis_vals,D_Date_Mar_1);plot(X_axis_vals,A_Date_Mar_1); % case 5
%% case 6: launch 15 march 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Mar_15(i,j),arrival_hyperbolic_excess_velocity_Date_Mar_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[F],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Mar_15 = departure_hyperbolic_excess_velocity_Date_Mar_15;D_Date_Mar_15 = Departure_Date_Mar_15';D_Date_Mar_15 = (D_Date_Mar_15(:))';
Arrival_Date_Mar_15 = arrival_hyperbolic_excess_velocity_Date_Mar_15;A_Date_Mar_15 = Arrival_Date_Mar_15';A_Date_Mar_15 = (A_Date_Mar_15(:))';
% plot(X_axis_vals,D_Date_Mar_15);plot(X_axis_vals,A_Date_Mar_15); % case 6
%% case 7: launch 1 april 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Apr_1(i,j),arrival_hyperbolic_excess_velocity_Date_Apr_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[G],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Apr_1 = departure_hyperbolic_excess_velocity_Date_Apr_1;D_Date_Apr_1 = Departure_Date_Apr_1';D_Date_Apr_1 = (D_Date_Apr_1(:))';
Arrival_Date_Apr_1 = arrival_hyperbolic_excess_velocity_Date_Apr_1;A_Date_Apr_1 = Arrival_Date_Apr_1';A_Date_Apr_1 = (A_Date_Apr_1(:))';
% plot(X_axis_vals,D_Date_Apr_1);plot(X_axis_vals,A_Date_Apr_1); % case 7
%% case 8: launch 15 april 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Apr_15(i,j),arrival_hyperbolic_excess_velocity_Date_Apr_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[H],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Apr_15 = departure_hyperbolic_excess_velocity_Date_Apr_15;D_Date_Apr_15 = Departure_Date_Apr_15';D_Date_Apr_15 = (D_Date_Apr_15(:))';
Arrival_Date_Apr_15 = arrival_hyperbolic_excess_velocity_Date_Apr_15;A_Date_Apr_15 = Arrival_Date_Apr_15';A_Date_Apr_15 = (A_Date_Apr_15(:))';
% plot(X_axis_vals,D_Date_Apr_15);plot(X_axis_vals,A_Date_Apr_15); % case 8
%% case 9: launch 1 may 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_May_1(i,j),arrival_hyperbolic_excess_velocity_Date_May_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[I],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_May_1 = departure_hyperbolic_excess_velocity_Date_May_1;D_Date_May_1 = Departure_Date_May_1';D_Date_May_1 = (D_Date_May_1(:))';
Arrival_Date_May_1 = arrival_hyperbolic_excess_velocity_Date_May_1;A_Date_May_1 = Arrival_Date_May_1';A_Date_May_1 = (A_Date_May_1(:))';
% plot(X_axis_vals,D_Date_May_1);plot(X_axis_vals,A_Date_May_1); % case 9
%% case 10: launch 15 may 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_May_15(i,j),arrival_hyperbolic_excess_velocity_Date_May_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[J],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_May_15 = departure_hyperbolic_excess_velocity_Date_May_15;D_Date_May_15 = Departure_Date_May_15';D_Date_May_15 = (D_Date_May_15(:))';
Arrival_Date_May_15 = arrival_hyperbolic_excess_velocity_Date_May_15;A_Date_May_15 = Arrival_Date_May_15';A_Date_May_15 = (A_Date_May_15(:))';
% plot(X_axis_vals,D_Date_May_15);plot(X_axis_vals,A_Date_May_15); % case 10
%% case 11: launch 1 june 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_June_1(i,j),arrival_hyperbolic_excess_velocity_Date_June_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[K],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_June_1 = departure_hyperbolic_excess_velocity_Date_June_1;D_Date_June_1 = Departure_Date_June_1';D_Date_June_1 = (D_Date_June_1(:))';
Arrival_Date_June_1 = arrival_hyperbolic_excess_velocity_Date_June_1;A_Date_June_1 = Arrival_Date_June_1';A_Date_June_1 = (A_Date_June_1(:))';
%plot(X_axis_vals,D_Date_June_1);plot(X_axis_vals,A_Date_June_1); % case 11
%% case 12: launch 15 june 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_June_15(i,j),arrival_hyperbolic_excess_velocity_Date_June_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[L],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_June_15 = departure_hyperbolic_excess_velocity_Date_June_15;D_Date_June_15 = Departure_Date_June_15';D_Date_June_15 = (D_Date_June_15(:))';
Arrival_Date_June_15 = arrival_hyperbolic_excess_velocity_Date_June_15;A_Date_June_15 = Arrival_Date_June_15';A_Date_June_15 = (A_Date_June_15(:))';
% plot(X_axis_vals,D_Date_June_15);plot(X_axis_vals,A_Date_June_15); % case 12
%% case 13: launch 1 july 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_July_1(i,j),arrival_hyperbolic_excess_velocity_Date_July_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[M],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_July_1 = departure_hyperbolic_excess_velocity_Date_July_1;D_Date_July_1 = Departure_Date_July_1';D_Date_July_1 = (D_Date_July_1(:))';
Arrival_Date_July_1 = arrival_hyperbolic_excess_velocity_Date_July_1;A_Date_July_1 = Arrival_Date_July_1';A_Date_July_1 = (A_Date_July_1(:))';
% plot(X_axis_vals,D_Date_July_1);plot(X_axis_vals,A_Date_July_1); % case 13
%% case 14: launch 15 july 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_July_15(i,j),arrival_hyperbolic_excess_velocity_Date_July_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[N],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_July_15 = departure_hyperbolic_excess_velocity_Date_July_15;D_Date_July_15 = Departure_Date_July_15';D_Date_July_15 = (D_Date_July_15(:))';
Arrival_Date_July_15 = arrival_hyperbolic_excess_velocity_Date_July_15;A_Date_July_15 = Arrival_Date_July_15';A_Date_July_15 = (A_Date_July_15(:))';
% plot(X_axis_vals,D_Date_July_15);plot(X_axis_vals,A_Date_July_15); % case 14
%% case 15: launch 1 aug 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Aug_1(i,j),arrival_hyperbolic_excess_velocity_Date_Aug_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[O],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Aug_1 = departure_hyperbolic_excess_velocity_Date_Aug_1;D_Date_Aug_1 = Departure_Date_Aug_1';D_Date_Aug_1 = (D_Date_Aug_1(:))';
Arrival_Date_Aug_1 = arrival_hyperbolic_excess_velocity_Date_Aug_1;A_Date_Aug_1 = Arrival_Date_Aug_1';A_Date_Aug_1 = (A_Date_Aug_1(:))';
% plot(X_axis_vals,D_Date_Aug_1);plot(X_axis_vals,A_Date_Aug_1); % case 15
%% case 16: launch 15 aug 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Aug_15(i,j),arrival_hyperbolic_excess_velocity_Date_Aug_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[P],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Aug_15 = departure_hyperbolic_excess_velocity_Date_Aug_15;D_Date_Aug_15 = Departure_Date_Aug_15';D_Date_Aug_15 = (D_Date_Aug_15(:))';
Arrival_Date_Aug_15 = arrival_hyperbolic_excess_velocity_Date_Aug_15;A_Date_Aug_15 = Arrival_Date_Aug_15';A_Date_Aug_15 = (A_Date_Aug_15(:))';
%plot(X_axis_vals,D_Date_Aug_15);plot(X_axis_vals,A_Date_Aug_15); % case 16
%% case 17: launch 1 sept 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Sept_1(i,j),arrival_hyperbolic_excess_velocity_Date_Sept_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[Q],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Sept_1 = departure_hyperbolic_excess_velocity_Date_Sept_1;D_Date_Sept_1 = Departure_Date_Sept_1';D_Date_Sept_1 = (D_Date_Sept_1(:))';
Arrival_Date_Sept_1 = arrival_hyperbolic_excess_velocity_Date_Sept_1;A_Date_Sept_1 = Arrival_Date_Sept_1';A_Date_Sept_1 = (A_Date_Sept_1(:))';
% plot(X_axis_vals,D_Date_Sept_1);plot(X_axis_vals,A_Date_Sept_1); % case 17
%% case 18: launch 15 sept 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Sept_15(i,j),arrival_hyperbolic_excess_velocity_Date_Sept_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[R],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Sept_15 = departure_hyperbolic_excess_velocity_Date_Sept_15;D_Date_Sept_15 = Departure_Date_Sept_15';D_Date_Sept_15 = (D_Date_Sept_15(:))';
Arrival_Date_Sept_15 = arrival_hyperbolic_excess_velocity_Date_Sept_15;A_Date_Sept_15 = Arrival_Date_Sept_15';A_Date_Sept_15 = (A_Date_Sept_15(:))';
% plot(X_axis_vals,D_Date_Sept_15);plot(X_axis_vals,A_Date_Sept_15); % case 18
%% case 19: launch 1 oct 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Oct_1(i,j),arrival_hyperbolic_excess_velocity_Date_Oct_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[S],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Oct_1 = departure_hyperbolic_excess_velocity_Date_Oct_1;D_Date_Oct_1 = Departure_Date_Oct_1';D_Date_Oct_1 = (D_Date_Oct_1(:))';
Arrival_Date_Oct_1 = arrival_hyperbolic_excess_velocity_Date_Oct_1;A_Date_Oct_1 = Arrival_Date_Oct_1';A_Date_Oct_1 = (A_Date_Oct_1(:))';
% plot(X_axis_vals,D_Date_Oct_1);plot(X_axis_vals,A_Date_Oct_1); % case 19
%% case 20: launch 15 oct 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Oct_15(i,j),arrival_hyperbolic_excess_velocity_Date_Oct_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[T],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Oct_15 = departure_hyperbolic_excess_velocity_Date_Oct_15;D_Date_Oct_15 = Departure_Date_Oct_15';D_Date_Oct_15 = (D_Date_Oct_15(:))';
Arrival_Date_Oct_15 = arrival_hyperbolic_excess_velocity_Date_Oct_15;A_Date_Oct_15 = Arrival_Date_Oct_15';A_Date_Oct_15 = (A_Date_Oct_15(:))';
% plot(X_axis_vals,D_Date_Oct_15);plot(X_axis_vals,A_Date_Oct_15); % case 20
%% case 21: launch 1 nov 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Nov_1(i,j),arrival_hyperbolic_excess_velocity_Date_Nov_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[U],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Nov_1 = departure_hyperbolic_excess_velocity_Date_Nov_1;D_Date_Nov_1 = Departure_Date_Nov_1';D_Date_Nov_1 = (D_Date_Nov_1(:))';
Arrival_Date_Nov_1 = arrival_hyperbolic_excess_velocity_Date_Nov_1;A_Date_Nov_1 = Arrival_Date_Nov_1';A_Date_Nov_1 = (A_Date_Nov_1(:))';
% plot(X_axis_vals,D_Date_Nov_1);plot(X_axis_vals,A_Date_Nov_1); % case 21
%% case 22: launch 15 nov 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Nov_15(i,j),arrival_hyperbolic_excess_velocity_Date_Nov_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[V],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Nov_15 = departure_hyperbolic_excess_velocity_Date_Nov_15;D_Date_Nov_15 = Departure_Date_Nov_15';D_Date_Nov_15 = (D_Date_Nov_15(:))';
Arrival_Date_Nov_15 = arrival_hyperbolic_excess_velocity_Date_Nov_15;A_Date_Nov_15 = Arrival_Date_Nov_15';A_Date_Nov_15 = (A_Date_Nov_15(:))';
% plot(X_axis_vals,D_Date_Nov_15);plot(X_axis_vals,A_Date_Nov_15); % case 22
%% case 23: launch 1 dec 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Dec_1(i,j),arrival_hyperbolic_excess_velocity_Date_Dec_1(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[W],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Dec_1 = departure_hyperbolic_excess_velocity_Date_Dec_1;D_Date_Dec_1 = Departure_Date_Dec_1';D_Date_Dec_1 = (D_Date_Dec_1(:))';
Arrival_Date_Dec_1 = arrival_hyperbolic_excess_velocity_Date_Dec_1;A_Date_Dec_1 = Arrival_Date_Dec_1';A_Date_Dec_1 = (A_Date_Dec_1(:))';
% plot(X_axis_vals,D_Date_Dec_1);plot(X_axis_vals,A_Date_Dec_1); % case 23
%% case 24: launch 15 dec 2036
for i = 1:length(Years)
    for j = 1:length(Months)        
    [departure_hyperbolic_excess_velocity_Date_Dec_15(i,j),arrival_hyperbolic_excess_velocity_Date_Dec_15(i,j)]  = Interplanetary_Trajectory('Earth','Jupiter',[X],[Years(i) Months(j) 1 0 00 00]);
    end
end
Departure_Date_Dec_15 = departure_hyperbolic_excess_velocity_Date_Dec_15;D_Date_Dec_15 = Departure_Date_Dec_15';D_Date_Dec_15 = (D_Date_Dec_15(:))';
Arrival_Date_Dec_15 = arrival_hyperbolic_excess_velocity_Date_Dec_15;A_Date_Dec_15 = Arrival_Date_Dec_15';A_Date_Dec_15 = (A_Date_Dec_15(:))';
% plot(X_axis_vals,D_Date_Dec_15);plot(X_axis_vals,A_Date_Dec_15); % case 24
%%
%%
figure(1)
subplot(2,3,1)
plot(X_axis_vals,D_Date_Jan_1);hold on;grid on;plot(X_axis_vals,A_Date_Jan_1); % case 1
plot(X_axis_vals,D_Date_Jan_15);plot(X_axis_vals,A_Date_Jan_15); % case 2
plot(X_axis_vals,D_Date_Feb_1);plot(X_axis_vals,A_Date_Feb_1); % case 3
plot(X_axis_vals,D_Date_Feb_15);plot(X_axis_vals,A_Date_Feb_15); % case 4
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ');
title('Jan to Feb 2036 launch')
% legend('Jan 1 Departure V_i_n_f','Jan 1 Arrival V_i_n_f',...
%        'Feb 1 Departure V_i_n_f','Feb 1 Arrival V_i_n_f')
% legend('Jan 1 Departure Hyperbolic Excess Velocity V_i_n_f','Jan 1 Arrival Hyperbolic Excess Velocity V_i_n_f',...
%     'Jan 15 Departure Hyperbolic Excess Velocity V_i_n_f','Jan 15 Arrival Hyperbolic Excess Velocity V_i_n_f',...
%     'Feb 1 Departure Hyperbolic Excess Velocity V_i_n_f','Feb 1 Arrival Hyperbolic Excess Velocity V_i_n_f',...
%     'Feb 15 Departure Hyperbolic Excess Velocity V_i_n_f','Feb 15 Arrival Hyperbolic Excess Velocity V_i_n_f')
subplot(2,3,2)
plot(X_axis_vals,D_Date_Mar_1);hold on;grid on;plot(X_axis_vals,A_Date_Mar_1); % case 5
plot(X_axis_vals,D_Date_Mar_15);plot(X_axis_vals,A_Date_Mar_15); % case 6
plot(X_axis_vals,D_Date_Apr_1);plot(X_axis_vals,A_Date_Apr_1); % case 7
plot(X_axis_vals,D_Date_Apr_15);plot(X_axis_vals,A_Date_Apr_15); % case 8
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('March to April 2036 launch')
subplot(2,3,3)
plot(X_axis_vals,D_Date_May_1);hold on;grid on;plot(X_axis_vals,A_Date_May_1); % case 9
plot(X_axis_vals,D_Date_May_15);plot(X_axis_vals,A_Date_May_15); % case 10
plot(X_axis_vals,D_Date_June_1);plot(X_axis_vals,A_Date_June_1); % case 11
plot(X_axis_vals,D_Date_June_15);plot(X_axis_vals,A_Date_June_15); % case 12
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('May to June 2036 launch')
subplot(2,3,4)
plot(X_axis_vals,D_Date_July_1);hold on;grid on;plot(X_axis_vals,A_Date_July_1); % case 13
plot(X_axis_vals,D_Date_July_15);plot(X_axis_vals,A_Date_July_15); % case 14
plot(X_axis_vals,D_Date_Aug_1);plot(X_axis_vals,A_Date_Aug_1); % case 15
plot(X_axis_vals,D_Date_Aug_15);plot(X_axis_vals,A_Date_Aug_15); % case 16
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('July to Aug 2036 launch')
subplot(2,3,5)
plot(X_axis_vals,D_Date_Sept_1);hold on;grid on;plot(X_axis_vals,A_Date_Sept_1); % case 17
plot(X_axis_vals,D_Date_Sept_15);plot(X_axis_vals,A_Date_Sept_15); % case 18
plot(X_axis_vals,D_Date_Oct_1);plot(X_axis_vals,A_Date_Oct_1); % case 19
plot(X_axis_vals,D_Date_Oct_15);plot(X_axis_vals,A_Date_Oct_15); % case 20
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('Sept to Oct 2036 launch')
subplot(2,3,6)
plot(X_axis_vals,D_Date_Nov_1);hold on;grid on;plot(X_axis_vals,A_Date_Nov_1); % case 21
plot(X_axis_vals,D_Date_Nov_15);plot(X_axis_vals,A_Date_Nov_15); % case 22
plot(X_axis_vals,D_Date_Dec_1);plot(X_axis_vals,A_Date_Dec_1); % case 23
plot(X_axis_vals,D_Date_Dec_15);plot(X_axis_vals,A_Date_Dec_15); % case 24
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('Nov to Dec 2036 launch')
figure(2)
subplot(1,3,1)
plot(X_axis_vals,D_Date_July_1);hold on;grid on;plot(X_axis_vals,A_Date_July_1); % case 13
plot(X_axis_vals,D_Date_July_15);plot(X_axis_vals,A_Date_July_15); % case 14
plot(X_axis_vals,D_Date_Aug_1);plot(X_axis_vals,A_Date_Aug_1); % case 15
plot(X_axis_vals,D_Date_Aug_15);plot(X_axis_vals,A_Date_Aug_15); % case 16
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('July to Aug 2036 launch')
subplot(1,3,2)
plot(X_axis_vals,D_Date_Sept_1);hold on;grid on;plot(X_axis_vals,A_Date_Sept_1); % case 17
plot(X_axis_vals,D_Date_Sept_15);plot(X_axis_vals,A_Date_Sept_15); % case 18
plot(X_axis_vals,D_Date_Oct_1);plot(X_axis_vals,A_Date_Oct_1); % case 19
plot(X_axis_vals,D_Date_Oct_15);plot(X_axis_vals,A_Date_Oct_15); % case 20
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('Sept to Oct 2036 launch')
subplot(1,3,3)
plot(X_axis_vals,D_Date_Nov_1);hold on;grid on;plot(X_axis_vals,A_Date_Nov_1); % case 21
plot(X_axis_vals,D_Date_Nov_15);plot(X_axis_vals,A_Date_Nov_15); % case 22
plot(X_axis_vals,D_Date_Dec_1);plot(X_axis_vals,A_Date_Dec_1); % case 23
plot(X_axis_vals,D_Date_Dec_15);plot(X_axis_vals,A_Date_Dec_15); % case 24
xlabel('SOI arrival at target planet (EARTH YEARS)');
ylabel('Departure and Arrival V_i_n_f (km/sec) ')
title('Nov to Dec 2036 launch')
figure(3)
%subplot(2,1,1)
plot(X_axis_vals,D_Date_July_1);hold on;grid on;plot(X_axis_vals,A_Date_July_1); % case 13
%plot(X_axis_vals,D_Date_July_15);plot(X_axis_vals,A_Date_July_15); % case 14
plot(X_axis_vals,D_Date_Aug_1);plot(X_axis_vals,A_Date_Aug_1); % case 15
%plot(X_axis_vals,D_Date_Aug_15);plot(X_axis_vals,A_Date_Aug_15); % case 16
% xlabel('SOI arrival at target planet (EARTH YEARS)');
xlabel('SOI arrival at Jupiter (Earth Years)');
ylabel('Earth Departure and Jupiter Arrival V_i_n_f (km/sec) ')
title('July to August 2036 launch')
legend('July 1 launch Departure V_i_n_f','July 1 launch Arrival V_i_n_f',...
      'August 1 launch Departure V_i_n_f','August 1 launch Arrival V_i_n_f')
%subplot(2,1,2)
figure(4)
plot(X_axis_vals,D_Date_Sept_1);hold on;grid on;plot(X_axis_vals,A_Date_Sept_1); % case 17
%plot(X_axis_vals,D_Date_Sept_15);plot(X_axis_vals,A_Date_Sept_15); % case 18
plot(X_axis_vals,D_Date_Oct_1);plot(X_axis_vals,A_Date_Oct_1); % case 19
%plot(X_axis_vals,D_Date_Oct_15);plot(X_axis_vals,A_Date_Oct_15); % case 20
% xlabel('SOI arrival at target planet (EARTH YEARS)');
% ylabel('Departure and Arrival V_i_n_f (km/sec) ')
xlabel('SOI arrival at Jupiter (Earth Years)');
ylabel('Earth Departure and Jupiter Arrival V_i_n_f (km/sec) ');
title('September to October 2036 launch');
legend('September 1 launch Departure V_i_n_f','September 1 launch Arrival V_i_n_f',...
      'October 1 launch Departure V_i_n_f','October 1 launch Arrival V_i_n_f');
%% EARTH DEPARTURE TRAJECTORY TILL EARTH SOI
% Assuming 0 deg inclination of deaprture orbit and hyperbola
R_earth = 149.6*10^(6); % km Earth orbital Radius
R_earth_radius = 6378.137; % km Earth Equatorial Radius
mu_earth = 398600.4418; %km^3/s^2 earths gravitaitonal parameter
mu_sun = 1.32712440018*10^(11); % km^3/s^2 sun gravitational parameter
rp = R_earth_radius+200; % km, initial parking orbit radius from center of the earth
v_departure_inf = 9:0.001:12; % km/sec, range of deaparture hyperbolic excess velocities
v_departure_p = sqrt(((v_departure_inf).^2) + ((2*mu_earth))/rp); % km/sec, speed of spacecraft at periapsis of the departure hyperbola
v_c = sqrt(mu_earth/rp); % km/sec, speed of sapcecraft in circular parking orbit
delta_v = v_departure_p-v_c; % km/sec, delta v requirement
eccentricity = 1+((rp*((v_departure_inf).^2))/mu_earth); % eccentricity of deaprture hyperbola 
figure(5)
plot(v_departure_inf,v_departure_p,'Linewidth',2);hold on;grid on;
plot(v_departure_inf,delta_v,'Linewidth',2)
plot(v_departure_inf,eccentricity,'Linewidth',2)
C3 = v_departure_inf.^2; % km^2/sec^2
plot(v_departure_inf,C3,'Linewidth',2);
xlabel('Departure v_i_n_f [km/sec]');ylabel('Depature v_p [km/sec]/Delta v [km/sec]/Hyperbolic eccentricity/Characteristic Energy [km^2/sec^2]')
legend('Depature v_p','Delta v','Eccentricity of Departure Hyperbola','Characteristic Energy')
title('Study of Departure trajectories from Earth to Jupiter')
yticks([10,20,30,40,50,60,70,80,90,100,110,120,130,140,160])
%% Jupiter flyby analysis
mu_jupiter = 126.686534*10^(6);% km^3/sec^2 gravitational prameter of jupiter
R_jupiter_radius = 71492; % km Jupiter Equatorial Radius
%%







