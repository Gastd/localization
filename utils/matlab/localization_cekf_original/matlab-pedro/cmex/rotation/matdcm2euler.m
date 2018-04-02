function [rpy] = dcm2euler(R)

%%% pitch ~= pi/2
% % calgary phd
% % pitch = atan2(-R(3,1),sqrt(R(3,2)^2+R(3,3)^2));
% roll = atan2(R(3,2),R(3,3));
% yaw = atan2(R(2,1),R(1,1));
% % titterton
% pitch = asin(-R(3,1));

pitch = asin(-R(3,1));
roll = atan2(R(3,2)/cos(pitch),R(3,3)/cos(pitch));
yaw = atan2(R(2,1)/cos(pitch),R(1,1)/cos(pitch));

rpy = [roll,pitch,yaw]';

% 
% if nargout == 1
%     varargout(1) = {[roll,pitch,yaw]'};
% end
% 
% if nargout == 3
%     varargout(1) = {roll};
%     varargout(2) = {pitch};
%     varargout(3) = {yaw};
% end
