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
real_data.gps_validpmeasure = gps_validpmeasure;
real_data.gps_validvmeasure = gps_validvmeasure;
real_data.sonar_range = sonar_range;
real_data.sonar_validmeasure = sonar_validmeasure;
real_data.X0 = X0;
real_data.P0 = P0;
