clc;clear;
start_time = juliandate(2022,7,10);
end_time = juliandate(2022,8,10);
time = start_time:end_time;
position_moon = zeros(length(time),3);
position_earth = zeros(length(time),3);
for i = 1:length(time)
    position_moon(i,:) = planetEphemeris(time(i),'Earth','Moon','432t');
end
for i = 1:length(time)
    position_earth(i,:) = planetEphemeris(time(i),'Moon','Earth','432t');
end
X_pos = position_moon(:,1);
Y_pos = position_moon(:,2);
Z_pos = position_moon(:,3);
X_e = position_earth(:,1);
Y_e = position_earth(:,2);
Z_e = position_earth(:,3);
figure(1)
for i = 2:length(time)
    
    plot3(X_pos(i),Y_pos(i),Z_pos(i),'r');hold on
    plot3(X_pos(i-1),Y_pos(i-1),Z_pos(i-1),'r')
    pause(0.001)

end
% plot3(X_e,Y_e,Z_e)