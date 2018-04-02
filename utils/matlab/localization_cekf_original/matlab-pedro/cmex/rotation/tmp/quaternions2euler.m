function [varargout] = quaternions2euler(q)

% make quaternion unit length:
q = q / norm(q);

R = quaternions2dcm(q);
[roll,pitch,yaw] = dcm2euler(R);

if nargout == 1
    varargout(1) = {[roll,pitch,yaw]'};
end

if nargout == 3
    varargout(1) = {roll};
    varargout(2) = {pitch};
    varargout(3) = {yaw};
end
