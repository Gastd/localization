function  [pose] = estimator_kf_state2pose(pose,kf_structure,flagpropagateuncertainty)

[rpy] = quaternions2euler(kf_structure.X(1:4));
pose.roll = rpy(1);
pose.pitch = rpy(2);
pose.yaw = rpy(3);
pose.q0 = kf_structure.X(1);
pose.q1 = kf_structure.X(2);
pose.q2 = kf_structure.X(3);
pose.q3 = kf_structure.X(4);
pose.x = kf_structure.X(5);
pose.y = kf_structure.X(6);
pose.z = kf_structure.X(7);
pose.dx_dt = kf_structure.X(8);
pose.dy_dt = kf_structure.X(9);
pose.dz_dt = kf_structure.X(10);

if flagpropagateuncertainty == 0
    return
end

% uses UT for uncertainty propagation from quaternions state
[Ym,Py] = unscented_transform(kf_structure.X(1:4),kf_structure.P(1:4,1:4),'quaternions2euler',zeros(4,1),[1; 1; 1]);
pose.rollvariance = Py(1,1);
pose.pitchvariance = Py(2,2);
pose.yawvariance = Py(3,3);

pose.q0variance = kf_structure.P(1,1);
pose.q1variance = kf_structure.P(2,2);
pose.q2variance = kf_structure.P(3,3);
pose.q3variance = kf_structure.P(4,4);
pose.xvariance = kf_structure.P(5,5);
pose.yvariance = kf_structure.P(6,6);
pose.zvariance = kf_structure.P(7,7);
pose.dx_dtvariance = kf_structure.P(8,8);
pose.dy_dtvariance = kf_structure.P(9,9);
pose.dz_dtvariance = kf_structure.P(10,10);

