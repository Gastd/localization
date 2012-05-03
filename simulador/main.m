%%% Simulador de estimadores de posicionamento tri-dimensional de um
%%% veiculo por um Sistema de Navegacao Inercial (INS). Os calculos serao
%%% realizados utilizando quaternios.

%% Initialization

clear all;
% close all;
clc;

if ispc
    addpath('..\cmex\sensors');
    addpath('..\cmex\rotation');
    addpath('..\cmex\localization');
else
    addpath('../cmex/sensors');
    addpath('../cmex/rotation');
    addpath('../cmex/localization');
end

%% Configuration parameters

% Estimator choice
flagestimatorTRIAD = 0; %TRIAD (attitude) only estimator
flagestimatorRungeKutta = 0; %4th order Runge-Kutta position estimation
flagestimatorEKF2 = 0; % Extended Kalman Filter (attitude+position)
flagestimatorCEKF = 0; % Correlated measurement Extended Kalman Filter (best) (attitude+position)
flagestimatorEKF3 = 0; % EKF with TRIAD uncertinity modulation, and other goodies (testing)
flagestimatorEKF_Decoupled = 1; % EKF decoupled magnetometer/accelerometer measurements (testing)
flagestimatorUKF2 = 0; % Sigma-point UKF (attitude+position)

% Should the system estimate the accelerometer bias?
flagkf_estimateaccelerometerbias = 1;

% Should the system estimate the gyrometer bias?
flagkf_estimategyrometerbias = 0; % Not fully implemented yet!

% Should we plot the attitute in real time? (slow)
plotRealtimeData = 0;

% Cleanup at the end?
flagcleanup = 1;

% What should be plotted at the end?
flagplot3sigma = 0;
plotQuaternions = 1;
plotAttitude = 1;
plotPosition = 1;
plotSpeed = 1;
plotBiases = 1;
plotSensorData = 01;
plotPosition3D = 0;

% Simulation or real data?
flaguserealdata = 1; % 0 -> simulation, 1 -> real data

% In case of a simulation, should we use noise?
flagnoise = 1;

% Ground truth available, if real data??
flaggroundtruth = 1;

% If real data, specify the data file
if ispc
    realdatafilename = 'dados\anu_localization.mat'; %PC
else
    realdatafilename = 'dados/anu_localization.mat'; %UNIX / Mac
end

% Should we generate data before (with the original data) or generate on
% the fly?
flaggeneratedatabefore = 1;

% What trajectory should be used in case of simulation?
%trajectory_name = 'trajectory_figure8';
trajectory_name = 'trajectory_helix';
%trajectory_name = 'trajectory_hovering';

%% System initialization

if flaguserealdata
    data = read_data(realdatafilename);
    if flaggroundtruth
        real_pose = data.real_pose;
    end
else
    if flaggeneratedatabefore
        [data, real_pose] = generate_trajectory(trajectory_name, flagnoise);
    else
        [data] = simulate_trajectory(trajectory_name, flagnoise, 0);
    end
end

% Aux variables
Ts = data.Ts;
t = data.t;
G = data.local_gravity';
M = data.local_magnetic';

% Measurement data structures
if flaguserealdata
    flagnoise = 1;
end

measurements = struct(	'sonar',repmat(sonar(0, flagnoise),length(t),1),...
                        'imu',repmat(imu(0, G, flagnoise),length(t),1),...
                        'gps',repmat(gps(0, flagnoise),length(t),1),...
                        'magnetometer',repmat(magnetometer(0, M, flagnoise),length(t),1));

% For Runge-Kutta and TRIAD
initial_estimates = repmat(pose_create,length(t),1);
initial_estimates(1).x = data.X0(5);
initial_estimates(1).y = data.X0(6);
initial_estimates(1).z = data.X0(7);
initial_estimates(1).dx_dt = data.X0(8);
initial_estimates(1).dy_dt = data.X0(9);
initial_estimates(1).dz_dt = data.X0(10);
initial_estimates(1).q0 = data.X0(1);
initial_estimates(1).q1 = data.X0(2);
initial_estimates(1).q2 = data.X0(3);
initial_estimates(1).q3 = data.X0(4);
[rpy] = quaternions2euler(data.X0(1:4));
initial_estimates(1).roll = rpy(1);
initial_estimates(1).pitch = rpy(2);
initial_estimates(1).yaw = rpy(3);
initial_estimates(1).rollvariance = (4*pi/180)^2;
initial_estimates(1).pitch = (4*pi/180)^2;
initial_estimates(1).yawvariance = (4*pi/180)^2;
initial_estimates(1).xvariance = (1)^2;
initial_estimates(1).yvariance = (1)^2;
initial_estimates(1).zvariance = (1)^2;
initial_estimates(1).dx_dtvariance = (0.1)^2;
initial_estimates(1).dy_dtvariance = (0.1)^2;
initial_estimates(1).dz_dtvariance = (0.1)^2;

% For visualization
if plotRealtimeData
    vehiclestate.x = data.X0(5);
    vehiclestate.y = data.X0(6);
    vehiclestate.z = data.X0(7);
    [rpy] = quaternions2euler([data.X0(1) data.X0(2) data.X0(3) data.X0(4)]);
    vehiclestate.roll = rpy(1);
    vehiclestate.pitch = rpy(2);
    vehiclestate.yaw = rpy(3);
    position_estimates.x = zeros(1, length(t));
    position_estimates.y = zeros(1, length(t));
    position_estimates.z = zeros(1, length(t));

    figure; plot3(vehiclestate.y,vehiclestate.x,-vehiclestate.z); hold on; grid on;
    hvehicle = vehicle_draw(vehiclestate);
    dx = 3000;
    dy = 3000;
    dz = 1000;
    axis([vehiclestate.x-dx vehiclestate.x+dx vehiclestate.y-dy vehiclestate.y+dy vehiclestate.z-dz vehiclestate.z+dz]);
    xlabel('x (east) [m]'); ylabel('y (north) [m]'); zlabel('z (up) [m]');
    set(gca,'DataAspectRatio',[1 1 1]);
end

%% Filter initialization

if flagestimatorEKF2
    [ekf2_structure] = localization_filter_init(flagkf_estimateaccelerometerbias, flagkf_estimategyrometerbias);
    ekf2_structure.X = data.X0;
    ekf2_structure.P = data.P0;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_ekf2 = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_ekf2 = repmat(parameter_estimates_ekf2,length(t),1);
        ekf2_structure.X(11) = 0; % Initial state estimate for X accelerometer bias
        ekf2_structure.X(12) = 0; % Initial state estimate for Y accelerometer bias
        ekf2_structure.X(13) = 0; % Initial state estimate for Z accelerometer bias
        ekf2_structure.P(11,11) = 0.0;
        ekf2_structure.P(12,12) = 0.0;
        ekf2_structure.P(13,13) = 0.0;
    end
    ekf2_structure.Preset = ekf2_structure.P;
    pose_estimates_ekf2 = repmat(pose_create,length(t),1);    
end

if flagestimatorEKF_Decoupled
    [ekf_decoupled_structure] = localization_filter_init(flagkf_estimateaccelerometerbias, flagkf_estimategyrometerbias);
    ekf_decoupled_structure.X = data.X0;
    ekf_decoupled_structure.P = data.P0;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_ekf_decoupled = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_ekf_decoupled = repmat(parameter_estimates_ekf_decoupled,length(t),1);
        ekf_decoupled_structure.X(11) = 0; % Initial state estimate for X accelerometer bias
        ekf_decoupled_structure.X(12) = 0; % Initial state estimate for Y accelerometer bias
        ekf_decoupled_structure.X(13) = 0; % Initial state estimate for Z accelerometer bias
        ekf_decoupled_structure.P(11,11) = 0.0;
        ekf_decoupled_structure.P(12,12) = 0.0;
        ekf_decoupled_structure.P(13,13) = 0.0;
    end
    ekf_decoupled_structure.Preset = ekf_decoupled_structure.P;
    pose_estimates_ekf_decoupled = repmat(pose_create,length(t),1);    
end

if flagestimatorCEKF
    [cekf_structure] = localization_filter_init(flagkf_estimateaccelerometerbias, flagkf_estimategyrometerbias);
    cekf_structure.X = data.X0;
    cekf_structure.P = data.P0;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_cekf = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_cekf = repmat(parameter_estimates_cekf,length(t),1);
        cekf_structure.X(11) = 0; % Initial state estimate for X accelerometer bias
        cekf_structure.X(12) = 0; % Initial state estimate for Y accelerometer bias
        cekf_structure.X(13) = 0; % Initial state estimate for Z accelerometer bias
        cekf_structure.P(11,11) = 0.0;
        cekf_structure.P(12,12) = 0.0;
        cekf_structure.P(13,13) = 0.0;
    end
    cekf_structure.Preset = cekf_structure.P;
    pose_estimates_cekf = repmat(pose_create,length(t),1);
end

if flagestimatorEKF3
    flagkf_estimategyrometerbias = 0;
    [ekf3_structure] = localization_filter_init(flagkf_estimateaccelerometerbias, flagkf_estimategyrometerbias);
    ekf3_structure.X = data.X0;
    ekf3_structure.P = data.P0;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_ekf3 = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        if flagkf_estimategyrometerbias
            ekf3_structure.X(11) = 0; % Initial state estimate for X accelerometer bias
            ekf3_structure.X(12) = 0; % Initial state estimate for Y accelerometer bias
            ekf3_structure.X(13) = 0; % Initial state estimate for Z accelerometer bias
            ekf3_structure.P(11,11) = 0.0;
            ekf3_structure.P(12,12) = 0.0;
            ekf3_structure.P(13,13) = 0.0;
            
            ekf3_structure.X(14) = 0; % Initial state estimate for X gyrometer bias
            ekf3_structure.X(15) = 0; % Initial state estimate for Y gyrometer bias
            ekf3_structure.X(16) = 0; % Initial state estimate for Z gyrometer bias
            ekf3_structure.P(14,14) = 0.0;
            ekf3_structure.P(15,15) = 0.0;
            ekf3_structure.P(16,16) = 0.0;
        else
            ekf3_structure.X(11) = 0; % Initial state estimate for X gyrometer bias
            ekf3_structure.X(12) = 0; % Initial state estimate for Y gyrometer bias
            ekf3_structure.X(13) = 0; % Initial state estimate for Z gyrometer bias
            ekf3_structure.P(11,11) = 0.0;
            ekf3_structure.P(12,12) = 0.0;
            ekf3_structure.P(13,13) = 0.0;
        end
    else
        if flagkf_estimategyrometerbias
            ekf3_structure.X(11) = 0; % Initial state estimate for X gyrometer bias
            ekf3_structure.X(12) = 0; % Initial state estimate for Y gyrometer bias
            ekf3_structure.X(13) = 0; % Initial state estimate for Z gyrometer bias
            ekf3_structure.P(11,11) = 0.0;
            ekf3_structure.P(12,12) = 0.0;
            ekf3_structure.P(13,13) = 0.0;
        end
    end
    parameter_estimates_ekf3 = repmat(parameter_estimates_ekf3,length(t),1);
    ekf3_structure.Preset = ekf3_structure.P;
    pose_estimates_ekf3 = repmat(pose_create,length(t),1);
end

if flagestimatorUKF2
    [ukf2_structure] = localization_filter_init(flagkf_estimateaccelerometerbias, flagkf_estimategyrometerbias);
    ukf2_structure.X = data.X0;
    ukf2_structure.P = data.P0;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_ukf2 = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_ukf2 = repmat(parameter_estimates_ukf2,length(t),1);
        ukf2_structure.X(11) = 0; % Initial state estimate for X accelerometer bias
        ukf2_structure.X(12) = 0; % Initial state estimate for Y accelerometer bias
        ukf2_structure.X(13) = 0; % Initial state estimate for Z accelerometer bias
        ukf2_structure.P(11,11) = 0.0;
        ukf2_structure.P(12,12) = 0.0;
        ukf2_structure.P(13,13) = 0.0;
    end
    ukf2_structure.Preset = ukf2_structure.P;
    pose_estimates_ukf2 = repmat(pose_create,length(t),1);
end

if flagestimatorRungeKutta
    pose_estimates_runge_kutta = initial_estimates;
end

if flagestimatorTRIAD
    pose_estimates_triad = initial_estimates;
    pose_estimates_triad(1).x = 0.0;
    pose_estimates_triad(1).y = 0.0;
    pose_estimates_triad(1).z = 0.0;
end

clear initial_estimates;

%% Main filter

for n=1:length(t),
    if rem(n,10)==0
        fprintf('*** t = %f\n',t(n));
    end
    
    if ((~flaggeneratedatabefore) && (~flaguserealdata))
        [data] = simulate_trajectory(trajectory_name, flagnoise, n);
        %%% IMU measurements:
        measurements.imu(n).wx = data.wx;
        measurements.imu(n).wy = data.wy;
        measurements.imu(n).wz = data.wz;
        measurements.imu(n).ax = data.ax;
        measurements.imu(n).ay = data.ay;
        measurements.imu(n).az = data.az;
        
        %%% Magnetometer measurements:
        measurements.magnetometer(n).mx = data.mx;
        measurements.magnetometer(n).my = data.my;
        measurements.magnetometer(n).mz = data.mz;
        measurements.magnetometer(n).flagvalidmeasure = data.magnetometer_validmeasure;
        
        %%% Sonar measurements:
        measurements.sonar(n).range = data.sonar_range;
        measurements.sonar(n).flagvalidmeasure = data.sonar_validmeasure;
        
        %%% GPS measurements:
        measurements.gps(n).p(1) = data.gps_x;
        measurements.gps(n).p(2) = data.gps_y;
        measurements.gps(n).p(3) = data.gps_z;
        measurements.gps(n).flagvalidpmeasure = data.gps_validpmeasure;
        measurements.gps(n).v(1) = data.gps_vx;
        measurements.gps(n).v(2) = data.gps_vy;
        measurements.gps(n).v(3) = data.gps_vz;
        measurements.gps(n).flagvalidvmeasure = data.gps_validvmeasure;
    else
        %%% IMU measurements:
        measurements.imu(n).wx = data.wx(n);
        measurements.imu(n).wy = data.wy(n);
        measurements.imu(n).wz = data.wz(n);
        measurements.imu(n).ax = data.ax(n);
        measurements.imu(n).ay = data.ay(n);
        measurements.imu(n).az = data.az(n);
        
        %%% Magnetometer measurements:
        measurements.magnetometer(n).mx = data.mx(n);
        measurements.magnetometer(n).my = data.my(n);
        measurements.magnetometer(n).mz = data.mz(n);
        measurements.magnetometer(n).flagvalidmeasure = data.magnetometer_validmeasure(n);
        
        %%% Sonar measurements:
        measurements.sonar(n).range = data.sonar_range(n);
        measurements.sonar(n).flagvalidmeasure = data.sonar_validmeasure(n);
        
        %%% GPS measurements:
        measurements.gps(n).p(1) = data.gps_x(n);
        measurements.gps(n).p(2) = data.gps_y(n);
        measurements.gps(n).p(3) = data.gps_z(n);
        measurements.gps(n).flagvalidpmeasure = data.gps_validpmeasure(n);
        measurements.gps(n).v(1) = data.gps_vx(n);
        measurements.gps(n).v(2) = data.gps_vy(n);
        measurements.gps(n).v(3) = data.gps_vz(n);
        measurements.gps(n).flagvalidvmeasure = data.gps_validvmeasure(n);
    end

    %%% TRIAD-based attitude estimation
    if flagestimatorTRIAD
        if measurements.magnetometer(n).flagvalidmeasure
            if n ~= 1
                q_previous = [pose_estimates_triad(n-1).q0 pose_estimates_triad(n-1).q1 pose_estimates_triad(n-1).q2 pose_estimates_triad(n-1).q3]';
            else
                q_previous = [pose_estimates_triad(1).q0 pose_estimates_triad(1).q1 pose_estimates_triad(1).q2 pose_estimates_triad(1).q3]';
            end
            q = mexlocalization('TRIAD',measurements.imu(n), measurements.magnetometer(n), M, G, q_previous);
            pose_estimates_triad(n).q0 = q(1);
            pose_estimates_triad(n).q1 = q(2);
            pose_estimates_triad(n).q2 = q(3);
            pose_estimates_triad(n).q3 = q(4);
            [rpy] = quaternions2euler(q);
            pose_estimates_triad(n).roll = rpy(1);
            pose_estimates_triad(n).pitch = rpy(2);
            pose_estimates_triad(n).yaw = rpy(3);
        end
        vehiclestate.x = 0.0;
        vehiclestate.y = 0.0;
        vehiclestate.z = 0.0;
        position_estimates.x(n) = vehiclestate.x;
        position_estimates.y(n) = vehiclestate.y;
        position_estimates.z(n) = vehiclestate.z;
        [rpy] = quaternions2euler([pose_estimates_triad(n).q0 pose_estimates_triad(n).q1 pose_estimates_triad(n).q2 pose_estimates_triad(n).q3]);
        vehiclestate.roll = rpy(1);
        vehiclestate.pitch = rpy(2);
        vehiclestate.yaw = rpy(3);
    end
    
    %%% Pose estimator using 4th order Runge-Kutta:
    if flagestimatorRungeKutta
        if n == 1
            pose_estimates_runge_kutta(n) = localization('RUNGEKUTTA',pose_estimates_runge_kutta(1), measurements.imu(n), G, Ts); 
        else
            pose_estimates_runge_kutta(n) = localization('RUNGEKUTTA',pose_estimates_runge_kutta(n-1), measurements.imu(n), G, Ts);
        end
        vehiclestate.x = pose_estimates_runge_kutta(n).x;
        vehiclestate.y = pose_estimates_runge_kutta(n).y;
        vehiclestate.z = pose_estimates_runge_kutta(n).z;
        position_estimates.x(n) = vehiclestate.x;
        position_estimates.y(n) = vehiclestate.y;
        position_estimates.z(n) = vehiclestate.z;
        [rpy] = quaternions2euler([pose_estimates_runge_kutta(n).q0 pose_estimates_runge_kutta(n).q1 pose_estimates_runge_kutta(n).q2 pose_estimates_runge_kutta(n).q3]);
        vehiclestate.roll = rpy(1);
        vehiclestate.pitch = rpy(2);
        vehiclestate.yaw = rpy(3);
    end

    %%% Extended Kalman Filter based estimation
    if flagestimatorEKF_Decoupled
        [ekf_decoupled_structure] = localization('FILTER_EKF_DECOUPLED_PREDICTION', ekf_decoupled_structure, measurements.imu(n), G, Ts);
        [ekf_decoupled_structure] = localization('FILTER_EKF_DECOUPLED_CORRECTION', ekf_decoupled_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, Ts);

        pose_estimates_ekf_decoupled(n) = localization('FILTER_STATE2POSE',pose_estimates_ekf_decoupled(n), ekf_decoupled_structure, 1);
        if flagkf_estimateaccelerometerbias
            parameter_estimates_ekf_decoupled(n).bias_ax = ekf_decoupled_structure.X(11);
            parameter_estimates_ekf_decoupled(n).bias_ay = ekf_decoupled_structure.X(12);
            parameter_estimates_ekf_decoupled(n).bias_az = ekf_decoupled_structure.X(13);
        end
        vehiclestate.x = pose_estimates_ekf_decoupled(n).x;
        vehiclestate.y = pose_estimates_ekf_decoupled(n).y;
        vehiclestate.z = pose_estimates_ekf_decoupled(n).z;
        position_estimates.x(n) = vehiclestate.x;
        position_estimates.y(n) = vehiclestate.y;
        position_estimates.z(n) = vehiclestate.z;
        vehiclestate.roll = pose_estimates_ekf_decoupled(n).roll;
        vehiclestate.pitch = pose_estimates_ekf_decoupled(n).pitch;
        vehiclestate.yaw = pose_estimates_ekf_decoupled(n).yaw;
    end
    
    %%% Extended Kalman Filter based estimation
    if flagestimatorEKF2
        [ekf2_structure] = mexlocalization('FILTER_EKF2_PREDICTION', ekf2_structure, measurements.imu(n), G, Ts);
        [ekf2_structure] = mexlocalization('FILTER_EKF2_CORRECTION', ekf2_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, Ts);

        pose_estimates_ekf2(n) = localization('FILTER_STATE2POSE',pose_estimates_ekf2(n), ekf2_structure, 1);
        if flagkf_estimateaccelerometerbias
            parameter_estimates_ekf2(n).bias_ax = ekf2_structure.X(11);
            parameter_estimates_ekf2(n).bias_ay = ekf2_structure.X(12);
            parameter_estimates_ekf2(n).bias_az = ekf2_structure.X(13);
        end
        vehiclestate.x = pose_estimates_ekf2(n).x;
        vehiclestate.y = pose_estimates_ekf2(n).y;
        vehiclestate.z = pose_estimates_ekf2(n).z;
        position_estimates.x(n) = vehiclestate.x;
        position_estimates.y(n) = vehiclestate.y;
        position_estimates.z(n) = vehiclestate.z;
        vehiclestate.roll = pose_estimates_ekf2(n).roll;
        vehiclestate.pitch = pose_estimates_ekf2(n).pitch;
        vehiclestate.yaw = pose_estimates_ekf2(n).yaw;
    end

    if flagestimatorCEKF
        [cekf_structure] = mexlocalization('FILTER_CEKF_PREDICTION', cekf_structure, measurements.imu(n), G, Ts);
        [cekf_structure] = mexlocalization('FILTER_CEKF_CORRECTION', cekf_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, Ts);

        pose_estimates_cekf(n) = localization('FILTER_STATE2POSE',pose_estimates_cekf(n), cekf_structure, 1);

        if flagkf_estimateaccelerometerbias
            parameter_estimates_cekf(n).bias_ax = cekf_structure.X(11);
            parameter_estimates_cekf(n).bias_ay = cekf_structure.X(12);
            parameter_estimates_cekf(n).bias_az = cekf_structure.X(13);
        end
        vehiclestate.x = pose_estimates_cekf(n).x;
        vehiclestate.y = pose_estimates_cekf(n).y;
        vehiclestate.z = pose_estimates_cekf(n).z;
        position_estimates.x(n) = vehiclestate.x;
        position_estimates.y(n) = vehiclestate.y;
        position_estimates.z(n) = vehiclestate.z;
        vehiclestate.roll = pose_estimates_cekf(n).roll;
        vehiclestate.pitch = pose_estimates_cekf(n).pitch;
        vehiclestate.yaw = pose_estimates_cekf(n).yaw;
    end
    
    if flagestimatorEKF3
        % Modulate the process covariance R based on the acceleration
%         accel_error = abs(norm(G) - sqrt(measurements.imu(n).ax.^2 + measurements.imu(n).ay.^2 + measurements.imu(n).az.^2));
%         ekf3_structure.R_convertedmeasurementtriad_ekf3 = exp(200*accel_error)*diag([0.01 0.01 0.01 0.01]);
%         if(ekf3_structure.R_convertedmeasurementtriad_ekf3(1,1) > 1e4)
%             ekf3_structure.R_convertedmeasurementtriad_ekf3 = 1e4*diag([0.01 0.01 0.01 0.01]);
%         end
        measurements.magnetometer(n).flagvalidmeasure = 0;
        ekf3_structure.R_convertedmeasurementtriad_ekf3 = 1e20*diag([0.01 0.01 0.01 0.01]);
        [ekf3_structure] = localization('FILTER_EKF3_PREDICTION', ekf3_structure, measurements.imu(n), G, Ts);
        [ekf3_structure] = localization('FILTER_EKF3_CORRECTION', ekf3_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, Ts);

        pose_estimates_ekf3(n) = localization('FILTER_STATE2POSE',pose_estimates_ekf3(n), ekf3_structure, 1);

        if flagkf_estimateaccelerometerbias
            parameter_estimates_ekf3(n).bias_ax = ekf3_structure.X(11);
            parameter_estimates_ekf3(n).bias_ay = ekf3_structure.X(12);
            parameter_estimates_ekf3(n).bias_az = ekf3_structure.X(13);
        end
        vehiclestate.x = pose_estimates_ekf3(n).x;
        vehiclestate.y = pose_estimates_ekf3(n).y;
        vehiclestate.z = pose_estimates_ekf3(n).z;
        position_estimates.x(n) = vehiclestate.x;
        position_estimates.y(n) = vehiclestate.y;
        position_estimates.z(n) = vehiclestate.z;
        vehiclestate.roll = pose_estimates_ekf3(n).roll;
        vehiclestate.pitch = pose_estimates_ekf3(n).pitch;
        vehiclestate.yaw = pose_estimates_ekf3(n).yaw;
    end
    
    if flagestimatorUKF2
        [ukf2_structure] = localization('FILTER_UKF2_PREDICTION', ukf2_structure, measurements.imu(n), G, Ts);
        [ukf2_structure] = localization('FILTER_UKF2_CORRECTION', ukf2_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, Ts);

        pose_estimates_ukf2(n) = localization('FILTER_STATE2POSE',pose_estimates_ukf2(n), ukf2_structure, 1);

        if flagkf_estimateaccelerometerbias
            parameter_estimates_ukf2(n).bias_ax = ukf2_structure.X(11);
            parameter_estimates_ukf2(n).bias_ay = ukf2_structure.X(12);
            parameter_estimates_ukf2(n).bias_az = ukf2_structure.X(13);
        end
        vehiclestate.x = pose_estimates_ukf2(n).x;
        vehiclestate.y = pose_estimates_ukf2(n).y;
        vehiclestate.z = pose_estimates_ukf2(n).z;
        position_estimates.x(n) = vehiclestate.x;
        position_estimates.y(n) = vehiclestate.y;
        position_estimates.z(n) = vehiclestate.z;
        vehiclestate.roll = pose_estimates_ukf2(n).roll;
        vehiclestate.pitch = pose_estimates_ukf2(n).pitch;
        vehiclestate.yaw = pose_estimates_ukf2(n).yaw;
    end
    
    if plotRealtimeData
        % Plot in ENU coordinates
        if rem(n,ceil(0.2/Ts))==2
            hvehicle = vehicle_draw(vehiclestate,hvehicle);
            if(n<100)
                plot3(position_estimates.y(1:n), position_estimates.x(1:n), -position_estimates.z(1:n));
            else
                plot3(position_estimates.y((n-100):n), position_estimates.x((n-100):n), -position_estimates.z((n-100):n));
            end
            drawnow; 
       end
    end
end

%% Data visualization

sigmacolor = [1 1 1]*0.80;
linewidth = 2;

% TRIAD Estimates
if flagestimatorTRIAD
    lineformat = 'b-';
    estimatorname = 'TRIAD';
    figure;
    I = find([measurements.magnetometer.flagvalidmeasure]~=0);
    subplot(311); hold on;
    data = [pose_estimates_triad.roll];
    plot(t(I),data(I)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('roll [deg]');
    title(['Attitude from: ',estimatorname]);
    subplot(312); hold on;
    data = [pose_estimates_triad.pitch];
    plot(t(I),data(I)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('pitch [deg]');
    subplot(313); hold on;
    data = [pose_estimates_triad.yaw];
    plot(t(I),data(I)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('yaw [deg]'); xlabel('t [s]');
    
    if flaguserealdata
        if flaggroundtruth
            subplot(311); hold on;
            e = [real_pose.roll]; e = atan2(sin(e),cos(e));
            plot(t,(e)*180/pi,'-.k','LineWidth',1);
            legend('Estimated data','Real data');
            subplot(312); hold on;
            e = [real_pose.pitch]; e = atan2(sin(e),cos(e));
            plot(t,(e)*180/pi,'-.k','LineWidth',1);
            legend('Estimated data','Real data');
            subplot(313); hold on;
            e = [real_pose.yaw]; e = atan2(sin(e),cos(e));
            plot(t,(e)*180/pi,'-.k','LineWidth',1);
            legend('Estimated data','Real data');
        end
    else if flaggeneratedatabefore %If we are simulating, we also plot the real information
            subplot(311); hold on;
            e = [real_pose.roll]; e = atan2(sin(e),cos(e));
            plot(t,(e)*180/pi,'-.k','LineWidth',1);
            legend('Estimated data','Real data');
            subplot(312); hold on;
            e = [real_pose.pitch]; e = atan2(sin(e),cos(e));
            plot(t,(e)*180/pi,'-.k','LineWidth',1);
            legend('Estimated data','Real data');
            subplot(313); hold on;
            e = [real_pose.yaw]; e = atan2(sin(e),cos(e));
            plot(t,(e)*180/pi,'-.k','LineWidth',1);
            legend('Estimated data','Real data');
        end
    end 
end

for nestimator=1:6
    switch nestimator
        case 1
            flagplot = flagestimatorEKF2;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_ekf2;
            if flagkf_estimateaccelerometerbias
                parameter_estimates = parameter_estimates_ekf2;
            end
            lineformat = 'b-';
            estimatorname = 'EKF2';
        case 2
            flagplot = flagestimatorCEKF;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_cekf;
            if flagkf_estimateaccelerometerbias
                parameter_estimates = parameter_estimates_cekf;
            end
            lineformat = 'r-';
            estimatorname = 'CEKF';
        case 3
            flagplot = flagestimatorRungeKutta;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_runge_kutta;
            lineformat = 'm-';
            estimatorname = 'Runge Kutta';
        case 4
            flagplot = flagestimatorUKF2;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_ukf2;
            if flagkf_estimateaccelerometerbias
                parameter_estimates = parameter_estimates_ukf2;
            end
            lineformat = 'g-';
            estimatorname = 'UKF2';            
        case 5
            flagplot = flagestimatorEKF3;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_ekf3;
            if flagkf_estimateaccelerometerbias
                parameter_estimates = parameter_estimates_ekf3;
            end
            lineformat = 'y-';
            estimatorname = 'EKF3';
        case 6
            flagplot = flagestimatorEKF_Decoupled;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_ekf_decoupled;
            if flagkf_estimateaccelerometerbias
                parameter_estimates = parameter_estimates_ekf_decoupled;
            end
            lineformat = 'm-';
            estimatorname = 'EKF Decoupled';
    end

    if plotQuaternions
        % Plot attitutes, in quaternions 
        figure;    
        subplot(411); hold on;
        title(['Attitude data (quaternions): ',estimatorname]);
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.q0variance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).q0 3*sigma(i-1)+pose_estimates(i-1).q0 -3*sigma(i-1)+pose_estimates(i-1).q0 -3*sigma(i)+pose_estimates(i).q0],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,([pose_estimates.q0]),lineformat,'LineWidth',linewidth);
        ylabel('_{q0}');
        xlabel('t [s]');

        subplot(412); hold on;
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.q1variance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).q1 3*sigma(i-1)+pose_estimates(i-1).q1 -3*sigma(i-1)+pose_estimates(i-1).q1 -3*sigma(i)+pose_estimates(i).q1],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,([pose_estimates.q1]),lineformat,'LineWidth',linewidth);
        ylabel('{q1}');
        xlabel('t [s]');

        subplot(413); hold on;
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.q2variance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).q2 3*sigma(i-1)+pose_estimates(i-1).q2 -3*sigma(i-1)+pose_estimates(i-1).q2 -3*sigma(i)+pose_estimates(i).q2],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,([pose_estimates.q2]),lineformat,'LineWidth',linewidth);
        ylabel('{q2}');
        xlabel('t [s]');

        subplot(414); hold on;
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.q3variance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).q3 3*sigma(i-1)+pose_estimates(i-1).q3 -3*sigma(i-1)+pose_estimates(i-1).q3 -3*sigma(i)+pose_estimates(i).q3],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,([pose_estimates.q3]),lineformat,'LineWidth',linewidth);
        ylabel('{q3}');
        xlabel('t [s]');
        
        if flaguserealdata
            if flaggroundtruth
                subplot(411); hold on;
                plot(t,([real_pose.q0]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(412); hold on;
                plot(t,([real_pose.q1]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(413); hold on;
                plot(t,([real_pose.q2]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(414); hold on;
                plot(t,([real_pose.q3]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                
            end
        else if flaggeneratedatabefore %If we are simulating, we also plot the real information
                subplot(411); hold on;
                plot(t,([real_pose.q0]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(412); hold on;
                plot(t,([real_pose.q1]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(413); hold on;
                plot(t,([real_pose.q2]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(414); hold on;
                plot(t,([real_pose.q3]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
            end
        end
    end
    
    if plotAttitude
        % Attitute graphs in degrees
        figure;
        title(['Attitude data (degrees): ',estimatorname]);
        subplot(311); hold on;
        e = [pose_estimates.roll]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,lineformat,'LineWidth',linewidth);    
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.rollvariance])*180/pi;
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+e(i)*180/pi 3*sigma(i-1)+e(i-1)*180/pi -3*sigma(i-1)+e(i-1)*180/pi -3*sigma(i)+e(i)*180/pi],sigmacolor,'EdgeColor',sigmacolor); end
        end
        ylabel('roll [deg]');

        subplot(312); hold on;
        e = [pose_estimates.pitch]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,lineformat,'LineWidth',linewidth);
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.pitchvariance])*180/pi;
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+e(i)*180/pi 3*sigma(i-1)+e(i-1)*180/pi -3*sigma(i-1)+e(i-1)*180/pi -3*sigma(i)+e(i)*180/pi],sigmacolor,'EdgeColor',sigmacolor); end
        end
        ylabel('pitch [deg]');

        subplot(313); hold on;
        e = [pose_estimates.yaw]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,lineformat,'LineWidth',linewidth);
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.yawvariance])*180/pi;
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+e(i)*180/pi 3*sigma(i-1)+e(i-1)*180/pi -3*sigma(i-1)+e(i-1)*180/pi -3*sigma(i)+e(i)*180/pi],sigmacolor,'EdgeColor',sigmacolor); end
        end
        ylabel('yaw [deg]'); xlabel('t [s]');

        if flaguserealdata
            if flaggroundtruth
                subplot(311); hold on;
                e = [real_pose.roll]; e = atan2(sin(e),cos(e));
                plot(t,(e)*180/pi,'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(312); hold on;
                e = [real_pose.pitch]; e = atan2(sin(e),cos(e));
                plot(t,(e)*180/pi,'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(313); hold on;
                e = [real_pose.yaw]; e = atan2(sin(e),cos(e));
                plot(t,(e)*180/pi,'-.k','LineWidth',1);
                legend('Estimated data','Real data');
            end
        else if flaggeneratedatabefore %If we are simulating, we also plot the real information
                subplot(311); hold on;
                e = [real_pose.roll]; e = atan2(sin(e),cos(e));
                plot(t,(e)*180/pi,'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(312); hold on;
                e = [real_pose.pitch]; e = atan2(sin(e),cos(e));
                plot(t,(e)*180/pi,'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(313); hold on;
                e = [real_pose.yaw]; e = atan2(sin(e),cos(e));
                plot(t,(e)*180/pi,'-.k','LineWidth',1);
                legend('Estimated data','Real data');
            end
        end
    end

    if plotPosition
        % Position data in meters
        figure;
        subplot(311); hold on;
        title(['Position data: ',estimatorname]);
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.xvariance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).x 3*sigma(i-1)+pose_estimates(i-1).x -3*sigma(i-1)+pose_estimates(i-1).x -3*sigma(i)+pose_estimates(i).x],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,[pose_estimates.x],lineformat,'LineWidth',linewidth);
        ylabel('x [m]');

        subplot(312); hold on;
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.yvariance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).y 3*sigma(i-1)+pose_estimates(i-1).y -3*sigma(i-1)+pose_estimates(i-1).y -3*sigma(i)+pose_estimates(i).y],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,[pose_estimates.y],lineformat,'LineWidth',linewidth);
        ylabel('y [m]');

        subplot(313); hold on;
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.zvariance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).z 3*sigma(i-1)+pose_estimates(i-1).z -3*sigma(i-1)+pose_estimates(i-1).z -3*sigma(i)+pose_estimates(i).z],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,[pose_estimates.z],lineformat,'LineWidth',linewidth);
        ylabel('z [m]'); xlabel('t [s]');
        
        if flaguserealdata
            if flaggroundtruth
                subplot(311); hold on;
                plot(t,([real_pose.x]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(312); hold on;
                plot(t,([real_pose.y]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(313); hold on;
                plot(t,([real_pose.z]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
            end
        else if flaggeneratedatabefore %If we are simulating, we also plot the real information
                subplot(311); hold on;
                plot(t,([real_pose.x]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(312); hold on;
                plot(t,([real_pose.y]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(313); hold on;
                plot(t,([real_pose.z]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
            end
        end
    end
    
    if plotSpeed
        % Speed data in meters/second
        figure;
        subplot(311); hold on;
        title(['Velocity data: ',estimatorname]);
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.dx_dtvariance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).dx_dt 3*sigma(i-1)+pose_estimates(i-1).dx_dt -3*sigma(i-1)+pose_estimates(i-1).dx_dt -3*sigma(i)+pose_estimates(i).dx_dt],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,[pose_estimates.dx_dt],lineformat,'LineWidth',linewidth);
        ylabel('v_x [m/s]'); % title(['Error in velocity: ',estimatorname]);

        subplot(312); hold on;
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.dy_dtvariance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).dy_dt 3*sigma(i-1)+pose_estimates(i-1).dy_dt -3*sigma(i-1)+pose_estimates(i-1).dy_dt -3*sigma(i)+pose_estimates(i).dy_dt],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,[pose_estimates.dy_dt],lineformat,'LineWidth',linewidth);
        ylabel('v_y [m/s]');

        subplot(313); hold on;
        if(flagplot3sigma)
            sigma = sqrt([pose_estimates.dz_dtvariance]);
            for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i)+pose_estimates(i).dz_dt 3*sigma(i-1)+pose_estimates(i-1).dz_dt -3*sigma(i-1)+pose_estimates(i-1).dz_dt -3*sigma(i)+pose_estimates(i).dz_dt],sigmacolor,'EdgeColor',sigmacolor); end
        end
        plot(t,[pose_estimates.dz_dt],lineformat,'LineWidth',linewidth);
        ylabel('v_z [m/s]'); xlabel('t [s]');
        
        if flaguserealdata
            if flaggroundtruth
                subplot(311); hold on;
                plot(t,([real_pose.dx_dt]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(312); hold on;
                plot(t,([real_pose.dy_dt]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(313); hold on;
                plot(t,([real_pose.dz_dt]),'-.k','LineWidth',1);
                legend('Estimated data','Real data'); 
            end
        else if flaggeneratedatabefore %If we are simulating, we also plot the real information
                subplot(311); hold on;
                plot(t,([real_pose.dx_dt]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(312); hold on;
                plot(t,([real_pose.dy_dt]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
                subplot(313); hold on;
                plot(t,([real_pose.dz_dt]),'-.k','LineWidth',1);
                legend('Estimated data','Real data');
            end
        end
    end
    
    if plotBiases
        if flagkf_estimateaccelerometerbias && (flagestimatorCEKF || flagestimatorEKF2)
            figure;
            subplot(311); hold on;
            title(['Accelerometer bias estimates: ',estimatorname]);
            plot(t,[parameter_estimates.bias_ax],lineformat,'LineWidth',linewidth);
            ylabel('b_a_x [m/s^2]');
            subplot(312); hold on;
            plot(t,[parameter_estimates.bias_ay],lineformat,'LineWidth',linewidth);
            ylabel('b_a_y [m/s^2]');
            subplot(313); hold on;
            plot(t,[parameter_estimates.bias_az],lineformat,'LineWidth',linewidth);
            ylabel('b_a_z [m/s^2]');
            xlabel('t [s]');
        end
    end
end

if plotSensorData
    % IMU Measurements
    figure;
    subplot(611); plot(t,[measurements.imu.ax],'b'); ylabel('a_x [m/s^2]');
    title('IMU measurements');
    subplot(612); plot(t,[measurements.imu.ay],'b'); ylabel('a_y [m/s^2]');
    subplot(613); plot(t,[measurements.imu.az],'b'); ylabel('a_z [m/s^2]');
    subplot(614); plot(t,[measurements.imu.wx]*180/pi,'b'); ylabel('w_x [deg/s]');
    subplot(615); plot(t,[measurements.imu.wy]*180/pi,'b'); ylabel('w_y [deg/s]');
    subplot(616); plot(t,[measurements.imu.wz]*180/pi,'b'); ylabel('w_z [deg/s]');
    xlabel('t [s]');

    % Magnetometer Measumerment
    figure;
    I = find([measurements.magnetometer.flagvalidmeasure]'~=0);
    plot_data = [measurements.magnetometer.mx];
    subplot(311); plot(t(I),plot_data(I),'b'); ylabel('m_x [T]');
    title('Magnetometer measurements');
    plot_data = [measurements.magnetometer.my];
    subplot(312); plot(t(I),plot_data(I),'b'); ylabel('m_y [T]');
    plot_data = [measurements.magnetometer.mz];
    subplot(313); plot(t(I),plot_data(I),'b'); ylabel('m_z [T]');
    xlabel('t [s]');

    % Sonar measurements
    figure;
    I = find([measurements.sonar.flagvalidmeasure]'~=0);
    plot_data = [measurements.sonar.range];
    plot(t(I),plot_data(1,I),'b'); ylabel('r [m]');
    title('Sonar measurements');

    % GPS Position measurements
    figure;
    I = find([measurements.gps.flagvalidpmeasure]'~=0);
    plot_data = [measurements.gps.p];
    subplot(311); plot(t(I),plot_data(1,I),'b.'); ylabel('x (north) [m]');
    title('GPS measurements: position');
    subplot(312); plot(t(I),plot_data(2,I),'b.'); ylabel('y (east) [m]');
    subplot(313); plot(t(I),plot_data(3,I),'b.'); ylabel('z (down) [m]');
    xlabel('t [s]');

    % GPS Velocity measurements
    figure;
    I = find([measurements.gps.flagvalidvmeasure]~=0);
    plot_data = [measurements.gps.v];
    subplot(311); plot(t(I),plot_data(1,I),'b.'); ylabel('v_x (north) [m/s]');
    title('GPS measurements: velocity');
    subplot(312); plot(t(I),plot_data(2,I),'b.'); ylabel('v_y (east) [m/s]');
    subplot(313); plot(t(I),plot_data(3,I),'b.'); ylabel('v_z (down) [m/s]');
    xlabel('t [s]');
end

% 3D plot of position estimates and GPS readings
% Note: we plot in the ENU system, so Z is up!
if plotPosition3D
    figure;
    title('3D Position data');
    Igps = find([measurements.gps.flagvalidpmeasure]~=0); plot_data = [measurements.gps.p];
    plot3(plot_data(2,Igps),plot_data(1,Igps),-1.0*plot_data(3,Igps),'k.'); hold on;
    if flagestimatorEKF2
        plot3([pose_estimates_ekf2.y],[pose_estimates_ekf2.x],-[pose_estimates_ekf2.z],'b');
    end
    if flagestimatorCEKF
        plot3([pose_estimates_cekf.y],[pose_estimates_cekf.x],-[pose_estimates_cekf.z],'r');
    end
    ylabel('y (north) [m]'); xlabel('x (east) [m]'); zlabel('z (up) [m]');
    set(gca,'DataAspectRatio',[1 1 1]);
end

%% Cleanup

if flagcleanup
    clear flagestimatorEKF3 flaggroundtruth I dx dy dz hvehicle plotRealtimeData vehiclestate flagcleanup position_estimates plotSensorData G e plot_data Igps M Ts cekf_structure ekf2_structure estimatorname flagestimatorCEKF flagestimatorEKF2 flagestimatorRungeKutta flagestimatorTRIAD flagkf_estimateaccelerometerbias flagnoise flagplot flagplot3sigma flaguserealdata lineformat linewidth n nestimator parameter_estimates plotAttitude plotBiases plotPosition plotPosition3D plotQuaternions plotSpeed pose_estimates realdatafilename rpy sigmacolor t trajectory_name flaggeneratedatabefore flagestimatorUKF2
end

return;
