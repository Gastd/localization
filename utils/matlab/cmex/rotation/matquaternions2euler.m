%Observar que na linha 1 coluna 2 h� um erro na Matriz do Padilha
function [rpy] = quaternions2euler(q)


% make quaternion unit length:
q = q / norm(q);

R = quaternions2dcm(q);
[rpy] = dcm2euler(R);
