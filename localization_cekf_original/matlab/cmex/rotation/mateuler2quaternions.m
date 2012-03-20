%Observar que na linha 1 coluna 2 há um erro na Matriz do Padilha
function q = euler2quaternions(rpy)

R = euler2dcm(rpy);
q = dcm2quaternions(R);
