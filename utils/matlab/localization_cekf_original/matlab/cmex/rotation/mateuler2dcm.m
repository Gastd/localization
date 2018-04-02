function [R] = euler2dcm(rpy)

roll = rpy(1);
pitch = rpy(2);
yaw = rpy(3);

R = zeros(3);
R(1,:) = [cos(pitch)*cos(yaw) -cos(roll)*sin(yaw)+sin(roll)*sin(pitch)*cos(yaw)  sin(roll)*sin(yaw)+cos(roll)*sin(pitch)*cos(yaw)];
R(2,:) = [cos(pitch)*sin(yaw)  cos(roll)*cos(yaw)+sin(roll)*sin(pitch)*sin(yaw) -sin(roll)*cos(yaw)+cos(roll)*sin(pitch)*sin(yaw)];
R(3,:) = [-sin(pitch)          sin(roll)*cos(pitch)                              cos(roll)*cos(pitch)];

