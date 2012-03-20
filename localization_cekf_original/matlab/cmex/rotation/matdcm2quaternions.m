function [q] = dcm2quaternions(R)

% improved code: Antonio Padilha
q0 = .5*sqrt(max(0,1+R(1,1)+R(2,2)+R(3,3)));
q1 = .5*sqrt(max(0,1+R(1,1)-R(2,2)-R(3,3)));
q2 = .5*sqrt(max(0,1-R(1,1)+R(2,2)-R(3,3)));
q3 = .5*sqrt(max(0,1-R(1,1)-R(2,2)+R(3,3)));

q1 = abs(q1)*sign(R(3,2)-R(2,3));
q2 = abs(q2)*sign(R(1,3)-R(3,1));
q3 = abs(q3)*sign(R(2,1)-R(1,2));

% traditional solution
% q0 = .5*sqrt(1 + R(1,1) + R(2,2) + R(3,3));
% q1 = (R(3,2) - R(2,3))/(4*q1);
% q2 = (R(1,3) - R(3,1))/(4*q1);
% q3 = (R(2,1) - R(1,2))/(4*q1);

q = [q0,q1,q2,q3]';