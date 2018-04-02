%função que calcula o quatérnio de rotação real acumulado, sem erros
function q = atitude_real(imu_data, q, T)

% cria a estrutura de distância angular
s = struct('x',0,'y',0,'z',0);

%preserva a atitude anterior
q_pre = q;

%atualiza a estrutura de distância ângular
s.x = imu_data.wx*T;
s.y = imu_data.wy*T;
s.z = imu_data.wz*T;

%Calcula o módulo
v = sqrt(s.x^2 + s.y^2 + s.z^2);


%Atualiza a atitude 
q(1) = cos(v/2)*q_pre(1) - sin(v/2)/v*(s.x*q_pre(2) + s.y*q_pre(3) + s.z*q_pre(4));

q(2) = cos(v/2)*q_pre(2) - sin(v/2)/v*(-s.x*q_pre(1) - s.z*q_pre(3) + s.y*q_pre(4));

q(3) = cos(v/2)*q_pre(3) - sin(v/2)/v*(-s.y*q_pre(1) + s.z*q_pre(2) - s.x*q_pre(4));

q(4) = cos(v/2)*q_pre(4) - sin(v/2)/v*(-s.z*q_pre(1) - s.y*q_pre(2) + s.x*q_pre(3));

