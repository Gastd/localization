%Observar que na linha 1 coluna 2 há um erro na Matriz do Padilha
function [rpy] = matquaternions2euler(q)


% make quaternion unit length:
q = q / norm(q);

R = quaternions2dcm(q);
[rpy] = dcm2euler(R);
