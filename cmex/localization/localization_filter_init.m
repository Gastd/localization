function  [kf_structure] = localization_filter_init(flagestimateaccelerometerbias, flagestimategyrometerbias)

if nargin == 0
    flagestimateaccelerometerbias = 0;
    flagestimategyrometerbias = 0;
end

if nargin == 1
    flagestimategyrometerbias = 0;
end

kf_structure.flagestimateaccelerometerbias = flagestimateaccelerometerbias;
kf_structure.flagestimategyrometerbias = flagestimategyrometerbias;

Nstates = 4 + 3 + 3;
if (flagestimateaccelerometerbias)
    Nstates = Nstates + 3;
end
if (flagestimategyrometerbias)
    Nstates = Nstates + 3;
end
kf_structure.Nstates = Nstates;

%Matrizes de covari�ncia
% Covariancia associada ao estado previsto
kf_structure.X      = zeros(Nstates,1);
kf_structure.P      = zeros(Nstates); %should be set by user.
kf_structure.Preset = zeros(Nstates); %should be set by user.
kf_structure.Xsigma = zeros(Nstates,2*Nstates+1);
kf_structure.Wsigma = zeros(1,2*Nstates+1);


% Covariancia associada ao processo. � uma matriz quadrada cujo n�mero de
% linhas � igual ao n�mero de linhas do processo.
kf_structure.Q_ekf = zeros(Nstates);
kf_structure.Q_ekf(1:4,1:4) = eye(4)*(0.01^2);
kf_structure.Q_ekf(5:7,5:7) = eye(3)*(0.1^2);
kf_structure.Q_ekf(8:10,8:10) = eye(3)*(0.01^2);
if (flagestimateaccelerometerbias)
    kf_structure.Q_ekf(11:13,11:13) = eye(3)*(0.001^2);
    if(flagestimategyrometerbias)
        kf_structure.Q_ekf(14:16,14:16) = eye(3)*(0.001^2);
    end
else
    if(flagestimategyrometerbias)
        kf_structure.Q_ekf(11:13,11:13) = eye(3)*(0.001^2);
    end
end

kf_structure.Q_ekf(1:4,1:4)   = kf_structure.Q_ekf(1:4,1:4)  * 1e0;
kf_structure.Q_ekf(5:10,5:10) = kf_structure.Q_ekf(5:10,5:10) * 1e-2;

kf_structure.Q_ukf = kf_structure.Q_ekf;
kf_structure.Q_cekf = kf_structure.Q_ekf;
kf_structure.Q_ekf3 = kf_structure.Q_ekf;
kf_structure.Q_ekf_decoupled = kf_structure.Q_ekf;

kf_structure.R_pseudomeasurementnorm = (1e-5)^2;
kf_structure.R_convertedmeasurementtriad_ekf = diag([0.01 0.01 0.01 0.01]);
kf_structure.R_convertedmeasurementtriad_ukf = diag([0.01 0.01 0.01 0.01]);
kf_structure.R_convertedmeasurementtriad_cekf = diag([0.01 0.01 0.01 0.01]*1e-1);
%kf_structure.R_convertedmeasurementtriad_cekf = diag([0.01 0.01 0.01 0.01]*1e2);
kf_structure.R_convertedmeasurementtriad_ekf3 = diag([0.01 0.01 0.01 0.01]);

kf_structure.S_previous_cekf = zeros(Nstates,4);
kf_structure.R_previous_cekf = eye(4);
kf_structure.A_imu_P_imu_cekf = zeros(Nstates,6);
kf_structure.innovation_previous_cekf = zeros(4,1);

kf_structure.S_previous_ekf3 = zeros(Nstates,4);
kf_structure.R_previous_ekf3 = eye(4);
kf_structure.A_imu_P_imu_ekf3 = zeros(Nstates,6);
kf_structure.innovation_previous_ekf3 = zeros(4,1);

kf_structure.ekf_decoupled_magnetometer_R = 1e-1*eye(3);
kf_structure.ekf_decoupled_accelerometer_R = 2*eye(3);
