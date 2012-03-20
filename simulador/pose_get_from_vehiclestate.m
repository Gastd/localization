function pose = pose_get_from_vehiclestate(vehiclestate)

pose = pose_create;

pose.x = vehiclestate.x;
pose.y = vehiclestate.y;
pose.z = vehiclestate.z;
pose.dx_dt = vehiclestate.dx_dt;
pose.dy_dt = vehiclestate.dy_dt;
pose.dz_dt = vehiclestate.dz_dt;
pose.pitch = vehiclestate.pitch;
pose.yaw = vehiclestate.yaw;
pose.roll = vehiclestate.roll;
q = euler2quaternions([vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw]');
pose.q0 = q(1);
pose.q1 = q(2);
pose.q2 = q(3);
pose.q3 = q(4);
