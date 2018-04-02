function  [kf_structure] = estimator_kf_init(flagestimateaccelerometerbias)

kf_structure.flagestimateaccelerometerbias = flagestimateaccelerometerbias;

Nstates = 4 + 3 + 3;
if (flagestimateaccelerometerbias)
    Nstates = Nstates + 3;
end
kf_structure.Nstates = Nstates;

%Matrizes de covariância
% Covariancia associada ao estado previsto
kf_structure.X      = zeros(Nstates,1);
kf_structure.P      = zeros(Nstates); %should be set by user.
kf_structure.Preset = zeros(Nstates); %should be set by user.
kf_structure.Xsigma = zeros(Nstates,2*Nstates+1);
kf_structure.Wsigma = zeros(1,2*Nstates+1);


% Covariancia associada ao processo. É uma matriz quadrada cujo número de
% linhas é igual ao número de linhas do processo.
kf_structure.Q_ekf = zeros(Nstates);
kf_structure.Q_ekf(1:4,1:4) = eye(4)*(0.01^2);
kf_structure.Q_ekf(5:7,5:7) = eye(3)*(0.1^2);
kf_structure.Q_ekf(8:10,8:10) = eye(3)*(0.01^2);
if (flagestimateaccelerometerbias)
    kf_structure.Q_ekf(11:13,11:13) = eye(3)*(0.001^2);
end
kf_structure.Q_ekf(1:4,1:4)   = kf_structure.Q_ekf(1:4,1:4)  * 1e0;
kf_structure.Q_ekf(5:10,5:10) = kf_structure.Q_ekf(5:10,5:10) * 1e-2;

kf_structure.Q_ukf = kf_structure.Q_ekf;
kf_structure.Q_cekf = kf_structure.Q_ekf;

kf_structure.R_pseudomeasurementnorm = (1e-5)^2;
kf_structure.R_convertedmeasurementtriad_ekf = diag([0.01 0.01 0.01 0.01]);
kf_structure.R_convertedmeasurementtriad_ukf = diag([0.01 0.01 0.01 0.01]);
kf_structure.R_convertedmeasurementtriad_cekf = diag([0.01 0.01 0.01 0.01]*1e-1);

kf_structure.S_previous_cekf = zeros(Nstates,4);
kf_structure.R_previous_cekf = eye(4);
kf_structure.A_imu_P_imu_cekf = zeros(Nstates,6);
kf_structure.innovation_previous_cekf = zeros(4,1);
