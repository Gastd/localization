function [real_data]  = read_data(filename)

%     The following data format is used in this program, keeping in mind
%     that they are ALL line vectors:
%     - Ts: sampling period, [s]
%     - t: time, [s]
%     - local_gravity: local gravity, in the format [acc_x acc_y acc_z] (NED coordinates), [m/s^2]
%     - local_magnetic: local magnetic field, in the format [mag_x mag_y mag_z] (NED coordinates), [uT]
%     - ax: acceleration (specific force) in the X axis of the body (e.g., aircraft), [m/s^2]
%     - ay: acceleration (specific force) in the Y axis of the body (e.g., aircraft), [m/s^2]
%     - az: acceleration (specific force) in the Z axis of the body (e.g., aircraft), [m/s^2]
%     - accelerometer_validmeasure: flag that sets if the measures from the accelerometer are valid (1 for valid, 0 for invalid).
%     - mx: magnetic field in the X axis of the body (e.g., aircraft), [uT]
%     - my: magnetic field in the Y axis of the body (e.g., aircraft), [uT]
%     - mz: magnetic field in the Z axis of the body (e.g., aircraft), [uT]
%     - magnetometer_validmeasure: flag that sets if the measures from the magnetometer are valid (1 for valid, 0 for invalid).
%     - wx: angular speed in the X axis of the body (e.g., aircraft), [rad/s]
%     - wy: angular speed in the Y axis of the body (e.g., aircraft), [rad/s]
%     - wz: angular speed in the Z axis of the body (e.g., aircraft), [rad/s]
%     - gyrometer_validmeasure: flag that sets if the measures from the gyrometers (angular speed) are valid (1 for valid, 0 for invalid).
%     - gps_x: X position from the GPS, in the LOCAL COORDINATE SYSTEM (NED coordinates, North). You get this by rotating the ECEF coordinates to the initial point. [m]
%     - gps_y: Y position from the GPS, in the LOCAL COORDINATE SYSTEM (NED coordinates, East). [m]
%     - gps_z: Z position from the GPS, in the LOCAL COORDINATE SYSTEM (NED coordinates, Down). [m]
%     - gps_validpmeasure: flag that sets if the position measurements from the GPS are valid (1 for valid, 0 for invalid).
%     - gps_vx: X velocity from the GPS, in the LOCAL COORDINATE SYSTEM (NED coordinates, North). You get this by rotating the velocities coordinates to the initial point. [m/s]
%     - gps_vy: Y velocity from the GPS, in the LOCAL COORDINATE SYSTEM (NED coordinates, East). [m/s]
%     - gps_vz: Z velocity from the GPS, in the LOCAL COORDINATE SYSTEM (NED coordinates, Down). [m/s]
%     - gps_validvmeasure: flag that sets if the velocity measurements from the GPS are valid (1 for valid, 0 for invalid).
%     - sonar_range: measurement from the sonar, if available [m]
%     - sonar_validmeasure: flag that sets if the sonar measurements are valid (1 for valid, 0 for invalid).
%     - altimeter: measurement from the altimeter, in height above ground [m]
%     - altimeter_validmeasure: flag that sets if the altimeter measurements are valid (1 for valid, 0 for invalid).
%     - X0: initial state estimate for the filter. 10x1 vector, with the
%     following structure: [quaternion_1 quaternion_2 quaternion_3
%     quaternion_4 position_x position_y position_z velocity_x velocity_y
%     velocity_z]. All in NED coordinates, with position in the local
%     coordinate system. Position in meters, velocity in meters/second.
%     - P0: initial error covariance matrix, [10x10], using the same structure as
%     the initial state estimate.

load(filename);

real_data.Ts = Ts;
real_data.t = t;
real_data.ax = ax;
real_data.ay = ay;
real_data.az = az;
real_data.accelerometer_validmeasure = accelerometer_validmeasure;
real_data.local_gravity = local_gravity;
real_data.wx = wx;
real_data.wy = wy;
real_data.wz = wz;
real_data.gyrometer_validmeasure = gyrometer_validmeasure;
real_data.mx = mx;
real_data.my = my;
real_data.mz = mz;
real_data.local_magnetic = local_magnetic;
real_data.magnetometer_validmeasure = magnetometer_validmeasure;
real_data.gps_x = gps_x;
real_data.gps_y = gps_y;
real_data.gps_z = gps_z;
real_data.gps_validpmeasure = gps_validpmeasure;
real_data.gps_vx = gps_vx;
real_data.gps_vy = gps_vy;
real_data.gps_vz = gps_vz;
real_data.gps_validvmeasure = gps_validvmeasure;
real_data.sonar_range = sonar_range;
real_data.sonar_validmeasure = sonar_validmeasure;

if(exist('real_data.X0', 'var'))
    real_data.X0 = X0;
else
    M = local_magnetic';
    G = local_gravity';

    % Perform a TRIAD estimate with the first data point
    mag = [mx(1, 1); my(1, 1); mz(1, 1)];

    %armazena a acelera��o (for�a espec�fica). A acerela��o � negada, pois leva-se em conta que os
    %acelerometros medem a for�a espec�fica sob sua massa e n�o a acelera��o
    a = zeros(3,1);
    a(1) = -ax(1, 1);
    a(2) = -ay(1, 1);
    a(3) = -az(1, 1);

    %normaliza��o dos vetores de acelera��o e de campo magn�tico
    mag = (mag/sqrt(mag(1)^2 + mag(2)^2 + mag(3)^2));
    a = a/sqrt(a(1)^2 + a(2)^2 + a(3)^2);
    G = G/sqrt(G(1)^2 + G(2)^2 + G(3)^2);
    M = M/sqrt(M(1)^2 + M(2)^2 + M(3)^2);


    % Inicializa as variaveis do corpo. sendo "a" as medidas do acelerometro e
    % "mag" as medidas do magnetometro. Lembrando que a fun��o cross(,)
    % representa produto vetorial entre dois vetores.
    aux = a + mag;
    I_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

    aux = cross(I_b,(a - mag));
    J_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

    K_b = cross(I_b,J_b);

    %Inicializa as variaveis do sistema inercial
    aux = G + M;
    I_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

    aux = cross(I_i,(G - M));
    J_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

    K_i = cross(I_i,J_i);

    % Calculo da matriz de rota��o pelo m�todo TRIAD melhorado
    R_i_b = ([I_b, J_b, K_b]*([I_i, J_i, K_i]'));

    %Cria��o do quat�rnio de rota��o. (Deve se observar que a matriz de
    %rota��o � do sistema de refer�ncia para o sistema do corpo a contr�rio do
    %que o Padilha escreveu. Dados para a confer�ncia est�o em
    %http://www.uel.br/proppg/semina/pdf/semina_28_1_22_19.pdf).
    R_b_i = R_i_b';
    q = dcm2quaternions(R_b_i);
    
    real_data.X0 = [q(1); q(2); q(3); q(4); gps_x(1); gps_y(1); gps_z(1); gps_vx(1); gps_vy(1); gps_vz(1);];
end

if(exist('real_data.P0', 'var'))
    real_data.P0 = P0;
else
    real_data.P0 = 1000.0*eye(10);
end

% Do we have ground truth?
if(exist('real_data_roll_radians', 'var'))
    real_pose = pose_create;
    real_pose.x = real_data_local_coordinate_system_north;
    real_pose.y = real_data_local_coordinate_system_east;
    real_pose.z = -1.0*real_data_local_coordinate_system_up;
    real_pose.dx_dt = real_data_northing_velocity_meters_second;
    real_pose.dy_dt = real_data_easting_velocity_meters_second;
    real_pose.dz_dt = real_data_down_velocity_meters_second;
    real_pose.pitch = real_data_pitch_radians;
    real_pose.yaw = real_data_yaw_radians;
    real_pose.roll = real_data_roll_radians;
    real_pose.q0 = zeros(length(real_data_roll_radians), 1);
    real_pose.q1 = zeros(length(real_data_roll_radians), 1);
    real_pose.q2 = zeros(length(real_data_roll_radians), 1);
    real_pose.q3 = zeros(length(real_data_roll_radians), 1);
    for i = 1:length(real_data_roll_radians)
        q = euler2quaternions([real_data_roll_radians(i), real_data_pitch_radians(i), real_data_yaw_radians(i)]');
        real_pose.q0(i) = q(1);
        real_pose.q1(i) = q(2);
        real_pose.q2(i) = q(3);
        real_pose.q3(i) = q(4);
    end
    real_data.real_pose = real_pose;
end
