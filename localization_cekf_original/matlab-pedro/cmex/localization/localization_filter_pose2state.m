%Converte e a pose do corpo para valores de suas variaveis de estado
%utilizadas durante a filtragem.

function  [kf_structure] = localization_filter_pose2state(pose,kf_structure,flagpropagateuncertainty)

kf_structure.X(1:4) = euler2quaternions([pose.roll, pose.pitch, pose.yaw]');
kf_structure.X(5) = pose.x;
kf_structure.X(6) = pose.y;
kf_structure.X(7) = pose.z;
kf_structure.X(8) = pose.dx_dt;
kf_structure.X(9) = pose.dy_dt;
kf_structure.X(10) = pose.dz_dt;

if flagpropagateuncertainty == 0
    return
end

% uses UT for uncertainty propagation to quaternions state
% in this case, covariances between q, p and v are lost.
[Ym,Py] = unscented_transform([pose.roll, pose.pitch, pose.yaw]',diag([pose.rollvariance, pose.pitchvariance, pose.yawvariance]),'euler2quaternions',[1; 1; 1],zeros(4,1));
kf_structure.P = zeros(length(kf_structure.X));
kf_structure.P(1:4,1:4) = Py;
kf_structure.P(5,5) = pose.xvariance;
kf_structure.P(6,6) = pose.yvariance;
kf_structure.P(7,7) = pose.zvariance;
kf_structure.P(8,8) = pose.dx_dtvariance;
kf_structure.P(9,9) = pose.dy_dtvariance;
kf_structure.P(10,10) = pose.dz_dtvariance;

