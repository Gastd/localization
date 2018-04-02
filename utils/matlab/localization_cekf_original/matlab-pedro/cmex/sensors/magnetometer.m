function magnetometermeasure = magnetometer(vehiclestate, M, flagnoise,mag_noise_multiplier)

magnetometermeasure = magnetometermeasure_init(flagnoise);

if isstruct(vehiclestate)==0
    return;
end

%Calcula a matriz de rota��o
R = euler2dcm([vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw]');

%Calculo do campo magnetico da terra em rela��o ao corpo a partir do
%refer�ncial inercial
bodymagfield = R'*M;
magnetometermeasure.mx = bodymagfield(1);
magnetometermeasure.my = bodymagfield(2);
magnetometermeasure.mz = bodymagfield(3);

%Os desvios padr�es dos magnet�metros s�o 0,0624 micro Tesla em x, 0,0274
%micro Tesla em y, 0,0546 micro Tesla em z
% Colocando ru�do
magnetometermeasure.mx = magnetometermeasure.mx + flagnoise*sqrt(magnetometermeasure.mxvariance*mag_noise_multiplier)*randn(1);
magnetometermeasure.my = magnetometermeasure.my + flagnoise*sqrt(magnetometermeasure.myvariance*mag_noise_multiplier)*randn(1);
magnetometermeasure.mz = magnetometermeasure.mz + flagnoise*sqrt(magnetometermeasure.mzvariance*mag_noise_multiplier)*randn(1);

magnetometermeasure.flagvalidmeasure = 1;
