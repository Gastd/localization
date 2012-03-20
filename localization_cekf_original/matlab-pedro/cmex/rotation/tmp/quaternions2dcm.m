function R = quaternions2dcm(q)

% make quaternion unit length:
q = q / norm(q);

%%% in acordance with Titterton & Weston, pg. 48, and Coutsias & Romero, last pg
R = 2*[(q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2)/2    (q(2)*q(3)-q(1)*q(4))                   (q(2)*q(4)+q(1)*q(3));
       (q(2)*q(3)+q(1)*q(4))                    (q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2)/2   (q(3)*q(4)-q(1)*q(2));
       (q(2)*q(4)-q(1)*q(3))                    (q(3)*q(4)+q(1)*q(2))                   (q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2)/2];

