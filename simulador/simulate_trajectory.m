function [simulated_trajectory_data_out] = simulate_trajectory(trajectory_name, flagnoise, n)

%% Simulation configuration parameters
sampling_frequency = 100; %simulation sampling frequency (Hz)
T_final = 20000.0; % total simulation time (seconds)

% Initial covariance multiplier
flaginitcovariance = 100;

% Simulated sensor subsampling frequencies
accelerometer_sampling_frequency = 100;
gyrometer_sampling_frequency = 100;
magnetometer_sampling_frequency = 100;
sonar_sampling_frequency = 0;
gps_position_sampling_frequency = 5;
gps_velocity_sampling_frequency = 5;

%Auxliary information
G = [0 0 9.781051614459725]'; % Local gravity
M = [19.9605 -7.6519 -9.8776]'; %Local magnetic field

%% Auxiliary variables and prealocation
persistent imu_measurement;
persistent magnetometer_measurement;
persistent sonar_measurement;
persistent gps_measurement;

Ts = 1/sampling_frequency;
t = (n-1)*Ts; % time
if n == 0
    t = 0;
end

%% If n=0, generate initial state estimate
if n == 0
    imu_measurement = imu(0, G, flagnoise);
    magnetometer_measurement = magnetometer(0, M, flagnoise);
    sonar_measurement = sonar(0, flagnoise);
    gps_measurement = gps(0, flagnoise);
    real_pose = repmat(pose_create, 1, 1);
    X0 = zeros(10, 1); %Initial state estimate. 4 quaternions + 3 positions + 3 velocities
    P0 = zeros(10, 10); %Initial error covariance matrix

    % Get simulated initial measurements for the vehicle
    [vehiclestate] = vehicle_getstate(t(1), trajectory_name);
    real_pose(1) = pose_get_from_vehiclestate(vehiclestate);
    imu_measurement = imu(vehiclestate, G, flagnoise);
    magnetometer_measurement = magnetometer(vehiclestate, M, flagnoise);
    gps_measurement = gps(vehiclestate, flagnoise);
    
    % Obtain the initial pose through a TRIAD estimate
    q_previous = [real_pose(1).q0 real_pose(1).q1 real_pose(1).q2 real_pose(1).q3]';
    q = mexlocalization('TRIAD',imu_measurement, magnetometer_measurement, M, G, q_previous);
    X0(1) = q(1); % Q0
    X0(2) = q(2); % Q1
    X0(3) = q(3); % Q2
    X0(4) = q(4); % Q3
    % And the position/speed through GPS
    X0(5) = gps_measurement.p(1); %X Position
    X0(6) = gps_measurement.p(2); %Y Position
    X0(7) = gps_measurement.p(3); %Z Position
    X0(8) = gps_measurement.v(1); %X Speed
    X0(9) = gps_measurement.v(2); %Y Speed
    X0(10) = gps_measurement.v(3); %Z Speed
    
    % Initial error covariance matrix
    rollvariance = flaginitcovariance*(4*pi/180)^2;
    pitchvariance = flaginitcovariance*(4*pi/180)^2;
    yawvariance = flaginitcovariance*(4*pi/180)^2;
    xvariance = flaginitcovariance*(1)^2;
    yvariance = flaginitcovariance*(1)^2;
    zvariance = flaginitcovariance*(1)^2;
    dx_dtvariance = flaginitcovariance*(0.1)^2;
    dy_dtvariance = flaginitcovariance*(0.1)^2;
    dz_dtvariance = flaginitcovariance*(0.1)^2;
    
    [~,Py] = unscented_transform([real_pose(1).roll, real_pose(1).pitch, real_pose(1).yaw]',diag([rollvariance, pitchvariance, yawvariance]),'euler2quaternions',[1; 1; 1],zeros(4,1));
    P0(1:4,1:4) = Py;
    P0(5,5) = xvariance;
    P0(6,6) = yvariance;
    P0(7,7) = zvariance;
    P0(8,8) = dx_dtvariance;
    P0(9,9) = dy_dtvariance;
    P0(10,10) = dz_dtvariance;
end

%% Measurment simulation
[vehiclestate] = vehicle_getstate(t, trajectory_name);

imu_measurement = imu(vehiclestate, G, flagnoise);
ax = imu_measurement.ax;
ay = imu_measurement.ay;
az = imu_measurement.az;
if (accelerometer_sampling_frequency ~= 0)
    if (mod(n-1, round(sampling_frequency/accelerometer_sampling_frequency)) == 0)
        accelerometer_validmeasure = 1;
    else
        accelerometer_validmeasure = 0;
    end
else
    accelerometer_validmeasure = 0;
end

wx = imu_measurement.wx;
wy = imu_measurement.wy;
wz = imu_measurement.wz;
if (gyrometer_sampling_frequency ~= 0)
    if (mod(n-1, round(sampling_frequency/gyrometer_sampling_frequency)) == 0)
        gyrometer_validmeasure = 1;
    else
        gyrometer_validmeasure = 0;
    end
else
    gyrometer_validmeasure = 0;
end

magnetometer_measurement = magnetometer(vehiclestate, M, flagnoise);
mx = magnetometer_measurement.mx;
my = magnetometer_measurement.my;
mz = magnetometer_measurement.mz;
if (magnetometer_sampling_frequency ~= 0)
    if (mod(n-1, round(sampling_frequency/magnetometer_sampling_frequency)) == 0)
        magnetometer_validmeasure = 1;
    else
        magnetometer_validmeasure = 0;
    end
else
    magnetometer_validmeasure = 0;
end

sonar_measurement = sonar(vehiclestate, flagnoise);
sonar_range = sonar_measurement.range;
if (sonar_sampling_frequency ~= 0)
    if (mod(n-1, round(sampling_frequency/sonar_sampling_frequency)) == 0)
        sonar_validmeasure = 1;
    else
        sonar_validmeasure = 0;
    end
else
    sonar_validmeasure = 0;
end

gps_measurement = gps(vehiclestate, flagnoise);    
gps_x = gps_measurement.p(1);
gps_y = gps_measurement.p(2);
gps_z = gps_measurement.p(3);
if (gps_position_sampling_frequency ~= 0)
    if (mod(n-1, round(sampling_frequency/gps_position_sampling_frequency)) == 0)
        gps_validpmeasure = 1;
    else
        gps_validpmeasure = 0;
    end
else
    gps_validpmeasure = 0;
end

gps_vx = gps_measurement.v(1);
gps_vy = gps_measurement.v(2);
gps_vz = gps_measurement.v(3);
if (gps_velocity_sampling_frequency ~= 0)
    if (mod(n-1, round(sampling_frequency/gps_velocity_sampling_frequency)) == 0)
        gps_validvmeasure = 1;
    else
        gps_validvmeasure = 0;
    end
else
    gps_validvmeasure = 0;
end

%% Output data
% Note: all of these are column vectors!
simulated_trajectory_data_out.t = t;
simulated_trajectory_data_out.ax = ax;
simulated_trajectory_data_out.ay = ay;
simulated_trajectory_data_out.az = az;
simulated_trajectory_data_out.accelerometer_validmeasure = accelerometer_validmeasure;
simulated_trajectory_data_out.wx = wx;
simulated_trajectory_data_out.wy = wy;
simulated_trajectory_data_out.wz = wz;
simulated_trajectory_data_out.gyrometer_validmeasure = gyrometer_validmeasure;
simulated_trajectory_data_out.mx = mx;
simulated_trajectory_data_out.my = my;
simulated_trajectory_data_out.mz = mz;
simulated_trajectory_data_out.magnetometer_validmeasure = magnetometer_validmeasure;
simulated_trajectory_data_out.gps_x = gps_x;
simulated_trajectory_data_out.gps_y = gps_y;
simulated_trajectory_data_out.gps_z = gps_z;
simulated_trajectory_data_out.gps_validpmeasure = gps_validpmeasure;
simulated_trajectory_data_out.gps_vx = gps_vx;
simulated_trajectory_data_out.gps_vy = gps_vy;
simulated_trajectory_data_out.gps_vz = gps_vz;
simulated_trajectory_data_out.gps_validvmeasure = gps_validvmeasure;
simulated_trajectory_data_out.sonar_range = sonar_range;
simulated_trajectory_data_out.sonar_validmeasure = sonar_validmeasure;

% If n=0, generate initial state estimate and the fixed data
if n == 0
    t = 0.0:Ts:T_final; % time
    simulated_trajectory_data_out.t = t;
    simulated_trajectory_data_out.Ts = Ts;
    simulated_trajectory_data_out.local_gravity = G';
    simulated_trajectory_data_out.local_magnetic = M';
    simulated_trajectory_data_out.X0 = X0;
    simulated_trajectory_data_out.P0 = P0;
end
