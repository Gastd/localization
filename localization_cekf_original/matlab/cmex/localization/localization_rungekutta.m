%Função que estima a posição, velocidade e atitude pelo método de Runge
%Kuta
function  [pose_estimator_state] = estimator_rungekutta(pose_estimator_state, imu_measurements, G, T)

%preserva a posição e suas 1ª  derivadas
V_pre = [pose_estimator_state.dx_dt; pose_estimator_state.dy_dt; pose_estimator_state.dz_dt];
P_pre = [pose_estimator_state.x; pose_estimator_state.y; pose_estimator_state.z];

%%% Calculo do quaternion
% ac garantirá que o modulo do quaternio será unitario
sigma_x = imu_measurements.wx*T;
sigma_y = imu_measurements.wy*T;
sigma_z = imu_measurements.wz*T;

sigma = sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2);

ac = 1 - (.5*sigma)^2/2 + (.5*sigma)^4/24;
as = .5*(1 - (.5*sigma)^2/2 + (.5*sigma)^4/24);

% Formatação do quatérnio
rk = [ac;
    as*sigma_x;
    as*sigma_y;
    as*sigma_z];

%% calcula a nova matriz de rotação
q_pre = euler2quaternions([pose_estimator_state.roll,pose_estimator_state.pitch,pose_estimator_state.yaw]');
q = multiplica_quaternion(q_pre,rk);
q = quaternions_correctsign(q, q_pre);

% now we check if we want q or -q
plus = q_pre - q;
minus = q_pre + q;
if (norm(plus) > norm(minus))
    q = - q;
end

%Calcula A no sistema inercial descontando a gravidade
A = quaternions2dcm(q_pre)*[imu_measurements.ax; imu_measurements.ay; imu_measurements.az] + G; % gravity compensation

%%% 4th order multivariable Runge-Kutta integration
% in acordance with http://www.myphysicslab.com/runge_kutta.html
V = V_pre + T*A;
k1 = T * ( A );
j1 = T * ( V );
k2 = T * ( A/2 );
j2 = T * ( V/2 );
k3 = T * ( A/2 );
j3 = T * ( V/2 );
k4 = T * ( A );
j4 = T * ( V );
V =  V_pre + k1/6 + k2/3 + k3/3 + k4/6;
P =  P_pre + j1/6 + j2/3 + j3/3 + j4/6;


%Armazena os dados calculados
[rpy] = quaternions2euler(q);
pose_estimator_state.roll = rpy(1);
pose_estimator_state.pitch = rpy(2);
pose_estimator_state.yaw = rpy(3);

pose_estimator_state.x = P(1);
pose_estimator_state.y = P(2);
pose_estimator_state.z = P(3);
pose_estimator_state.dx_dt = V(1);
pose_estimator_state.dy_dt = V(2);
pose_estimator_state.dz_dt = V(3);


