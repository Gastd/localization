function magnetometermeasure = magnetometro(vehiclestate, M, flagnoise)

magnetometermeasure = magnetometermeasure_init(flagnoise);

if isstruct(vehiclestate)==0
    return;
end

%Calcula a matriz de rotação
R = euler2dcm([vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw]');

%Calculo do campo magnetico da terra em relação ao corpo a partir do
%referêncial inercial
bodymagfield = R'*M;
magnetometermeasure.mx = bodymagfield(1);
magnetometermeasure.my = bodymagfield(2);
magnetometermeasure.mz = bodymagfield(3);

%Os desvios padrões dos magnetômetros são 0,0624 micro Tesla em x, 0,0274
%micro Tesla em y, 0,0546 micro Tesla em z
% Colocando ruído
magnetometermeasure.mx = magnetometermeasure.mx + flagnoise*sqrt(magnetometermeasure.mxvariance)*randn(1);
magnetometermeasure.my = magnetometermeasure.my + flagnoise*sqrt(magnetometermeasure.myvariance)*randn(1);
magnetometermeasure.mz = magnetometermeasure.mz + flagnoise*sqrt(magnetometermeasure.mzvariance)*randn(1);

magnetometermeasure.flagvalidmeasure = 1;
