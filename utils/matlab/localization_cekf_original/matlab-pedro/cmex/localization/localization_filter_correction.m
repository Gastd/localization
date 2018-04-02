%Função que a realiza a etapa de correção pelo Filtro de Kalman
function  [kf_structure] = localization_filter_correction(filtername, kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)

switch filtername
    case 'ukf2'
        [kf_structure] = localization_filter_correction_ukf2(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    case 'ekf2'
        [kf_structure] = localization_filter_correction_ekf2(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    case 'cekf'
        [kf_structure] = localization_filter_correction_cekf(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_correction_cekf(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)
%%% corrige usando medidas do algoritmo TRIAD: 
%%% Funçao de mediçao: [eye(4) zeros(4,6)]*X (a medida é dada pelo
%%% estimador TRIAD.
%%% Estimar a mediçao do quaternion assim como da matriz R usando UT;
if((magnetometermeasure.flagvalidmeasure))
    %%%
    X = kf_structure.X;
    P = kf_structure.P;
    if(kf_structure.flagestimateaccelerometerbias)
        H = [eye(4) zeros(4,9)];
    else
        H = [eye(4) zeros(4,6)];
    end
    [Ym,Py] = ConvertedMeasurementTRIAD('cekf',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
    v = (Ym - H*X);
    Rtilde = Py + kf_structure.R_convertedmeasurementtriad_cekf;
    Stilde = kf_structure.A_imu_P_imu_cekf * localization_filter_model_gtriad('dg_du_imu',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias,1)';
    K = P*(H')*inv(H*P*(H') + Rtilde);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
    kf_structure.R_previous_cekf = Rtilde;
    kf_structure.S_previous_cekf = Stilde;
    kf_structure.innovation_previous_cekf = (Ym - H*kf_structure.X);
end

%%% corrige usando medidas do GPS:
%%% Funçao de mediçao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
if gpsmeasure.flagvalidpmeasure
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
    else
        H = [zeros(3,4) eye(3,3) zeros(3,3)];
    end
    v = (gpsmeasure.p - H*X);
    K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do GPS:
%%% Funçao de mediçao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
if gpsmeasure.flagvalidvmeasure
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
    else
        H = [zeros(3,4) zeros(3,3) eye(3,3)];
    end
    v = (gpsmeasure.v - H*X);
    K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do sonar:
%%% Funçao de mediçao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Funçao de mediçao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
if sonarmeasure.flagvalidmeasure
    X = kf_structure.X;
    P = kf_structure.P;
    dHdX1 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(1);
    dHdX2 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(2);
    dHdX3 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(3);
    dHdX4 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(4);
    dHdX7 =  1/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
    if kf_structure.flagestimateaccelerometerbias
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0];
    else
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0];
    end
    v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
    K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige a norma do quaternion usando pseudo-mediçoes:
%%% Funçao de mediçao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Funçao de mediçao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
X = kf_structure.X;
P = kf_structure.P;
dHdX1 = 2*X(1);
dHdX2 = 2*X(2);
dHdX3 = 2*X(3);
dHdX4 = 2*X(4);
if kf_structure.flagestimateaccelerometerbias
    H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
else
    H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
end
v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
kf_structure.X = X + K*v;
kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_correction_ekf2(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)
%%% corrige usando medidas do algoritmo TRIAD: 
%%% Funçao de mediçao: [eye(4) zeros(4,6)]*X (a medida é dada pelo
%%% estimador TRIAD.
%%% Estimar a mediçao do quaternion assim como da matriz R usando UT;
if((magnetometermeasure.flagvalidmeasure == 1))
    %%% 
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        H = [eye(4) zeros(4,9)];
    else
        H = [eye(4) zeros(4,6)];
    end
    [Ym,Py] = ConvertedMeasurementTRIAD('ekf2',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
    v = (Ym - H*X);
    K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad_ekf);
    kf_structure.X = X + K*v;
    X = kf_structure.X;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%   Corrige a atitude utilizando diretamente as medidas do magnetometro,
%   sem passar pelo algoritmo TRIAD.
if((magnetometermeasure.flagvalidmeasure == 2)||(magnetometermeasure.flagvalidmeasure == 3))      
    %%%    
    X = kf_structure.X;
    P = kf_structure.P;
    
    %Calcula os valores numericos da matriz Jacobiana de medicao H e da
    %funcao de medicao h_est baseada na ultima estimativa de estado.
    
    %Modelo do magnetometro funcionando
    if magnetometermeasure.flagvalidmeasure == 2;    
        [H,h_est] = subs_EKF(kf_structure,M);    
    %Modelo do magnetometro com falha
    else
        H = zeros(3,13);  
        h_est = zeros(3,1);
    end
    
    %Vetor de medidas do magnetometro
    mag_meas = [magnetometermeasure.mx magnetometermeasure.my magnetometermeasure.mz]';
    %Matriz de covariancias da medida do magnetometro
    mag_cov = [magnetometermeasure.mxvariance 0 0; 0 magnetometermeasure.myvariance 0; 0 0 magnetometermeasure.mzvariance];
    
    %Termo de inovacao -> diferenca entre a medida do magnetometro e seu
    %valor esperado pela funcao h_mag.
    v = mag_meas - h_est;    
    %Ganho de Kalman
    K = P*(H')*inv(H*P*(H') + mag_cov);
    %Corrige a estimativa
    kf_structure.X = X + K*v;
    
    %X agora armazena a estimativa corrigida
    X = kf_structure.X;
    %Corrige o sinal do quaternio, pois q e -q representam a mesma
    %orientacao.
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    
    %Matriz de covariancias corrigida
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
    
end


%%% corrige usando medidas do GPS:
%%% Funçao de mediçao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
if gpsmeasure.flagvalidpmeasure
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
    else
        H = [zeros(3,4) eye(3,3) zeros(3,3)];
    end
    v = (gpsmeasure.p - H*X);
    K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do GPS:
%%% Funçao de mediçao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
if gpsmeasure.flagvalidvmeasure
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
    else
        H = [zeros(3,4) zeros(3,3) eye(3,3)];
    end
    v = (gpsmeasure.v - H*X);
    K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do sonar:
%%% Funçao de mediçao: z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Funçao de mediçao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
if(sonarmeasure.flagvalidmeasure)
    X = kf_structure.X;
    P = kf_structure.P;
    dHdX1 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(1);
    dHdX2 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(2);
    dHdX3 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(3);
    dHdX4 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(4);
    dHdX7 =  1/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
    if kf_structure.flagestimateaccelerometerbias
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0];
    else
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0];
    end
    v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
    K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige a norma do quaternion usando pseudo-mediçoes:
%%% Funçao de mediçao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Funçao de mediçao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
X = kf_structure.X;
P = kf_structure.P;
dHdX1 = 2*X(1);
dHdX2 = 2*X(2);
dHdX3 = 2*X(3);
dHdX4 = 2*X(4);
if kf_structure.flagestimateaccelerometerbias
    H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
else
    H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
end
v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
kf_structure.X = X + K*v;
kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_correction_ukf2(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)
%%% corrige usando medidas do algoritmo TRIAD: 
%%% Funçao de mediçao: [eye(4) zeros(4,6)]*X (a medida é dada pelo
%%% estimador TRIAD.
%%% Estimar a mediçao do quaternion assim como da matriz R usando UT;
if(isstruct(magnetometermeasure) & isstruct(imumeasure))
    %%% 
    X = kf_structure.X;
    P = kf_structure.P;
    H = [eye(4) zeros(4,6)];
    [Ym,Py] = ConvertedMeasurementTRIAD('ukf2',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
    v = (Ym - H*X);
    K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do sonar:
%%% Funçao de mediçao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Funçao de mediçao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
if(isstruct(sonarmeasure))
    X = kf_structure.X;
    P = kf_structure.P;
    dHdX1 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(1);
    dHdX2 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(2);
    dHdX3 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(3);
    dHdX4 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(4);
    dHdX7 =  1/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
    H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0];
    v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
    K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige a norma do quaternion usando pseudo-mediçoes:
%%% Funçao de mediçao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Funçao de mediçao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
X = kf_structure.X;
P = kf_structure.P;
dHdX1 = 2*X(1);
dHdX2 = 2*X(2);
dHdX3 = 2*X(3);
dHdX4 = 2*X(4);
H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
kf_structure.X = X + K*v;
kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ym,Py] = ConvertedMeasurementTRIAD(filtername, X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias)


switch filtername
    case 'ekf2'
        [Ym,Py] = ConvertedMeasurementTRIAD('cekf', X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias);
    case 'cekf'
        Ym =  localization_filter_model_gtriad('evaluate',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
        dg_du_imu = localization_filter_model_gtriad('dg_du_imu',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
        dg_du_mag = localization_filter_model_gtriad('dg_du_mag',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
        dg_dx = localization_filter_model_gtriad('dg_dx',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
        P_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
        P_mag = diag([magnetometermeasure.mxvariance; magnetometermeasure.myvariance; magnetometermeasure.mzvariance]);
        Py = dg_du_imu * P_imu * dg_du_imu' + dg_du_mag * P_mag * dg_du_mag' + dg_dx * Px * dg_dx';
    case 'ukf2'
        Xm = [magnetometermeasure.mx; magnetometermeasure.my; magnetometermeasure.mz; imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
        Pmag = diag([magnetometermeasure.mxvariance; magnetometermeasure.myvariance; magnetometermeasure.mzvariance]);
        Pimu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
        Px = [Pmag zeros(3,6); zeros(6,3) Pimu];
        DimX = 9;
        DimY = 4;
        kappa = 1;
        UtNSamples = 2*DimX+1;
        XSamples = zeros(DimX,UtNSamples);
        YSamples = zeros(DimY,UtNSamples);
        WSamples = zeros(1,UtNSamples);
        Q = sqrtm((DimX+kappa)*Px)';
        for i=1:DimX,
            WSamples(1,i)   = 1/(2*(DimX+kappa));
            WSamples(1,i+DimX) = 1/(2*(DimX+kappa));
            XSamples(:,i)   = Xm + Q(:,i);
            XSamples(:,i+DimX) = Xm - Q(:,i);
            X = XSamples(:,i);
            magnetometermeasure.mx = X(1); magnetometermeasure.my = X(2); magnetometermeasure.mz = X(3);
            imumeasure.ax = X(4); imumeasure.ay = X(5); imumeasure.az = X(6); imumeasure.wx = X(7); imumeasure.wy = X(8); imumeasure.wz = X(9);
            YSamples(:,i) =  localization_filter_model_gtriad('evaluate',imumeasure, magnetometermeasure, M, G, q_predicted);
            X = XSamples(:,i+DimX);
            magnetometermeasure.mx = X(1); magnetometermeasure.my = X(2); magnetometermeasure.mz = X(3);
            imumeasure.ax = X(4); imumeasure.ay = X(5); imumeasure.az = X(6); imumeasure.wx = X(7); imumeasure.wy = X(8); imumeasure.wz = X(9);
            YSamples(:,i+DimX) =  localization_filter_model_gtriad('evaluate',imumeasure, magnetometermeasure, M, G, q_predicted);
        end,
        WSamples(1,UtNSamples) = kappa/(DimX+kappa);
        XSamples(:,UtNSamples) = Xm;
        X = XSamples(:,UtNSamples);
        magnetometermeasure.mx = X(1); magnetometermeasure.my = X(2); magnetometermeasure.mz = X(3);
        imumeasure.ax = X(4); imumeasure.ay = X(5); imumeasure.az = X(6); imumeasure.wx = X(7); imumeasure.wy = X(8); imumeasure.wz = X(9);
        YSamples(:,UtNSamples) =  localization_filter_model_gtriad('evaluate',imumeasure, magnetometermeasure, M, G);
        Ym 	= zeros(DimY,1);
        Py  = zeros(DimY,DimY);
        for ns=1:UtNSamples
            Ym = Ym + WSamples(1,ns)*YSamples(:,ns);
        end,
        for ns=1:UtNSamples
            Py = Py + WSamples(1,ns)*(YSamples(:,ns)-Ym)*(YSamples(:,ns)-Ym)';
        end,
end


return;

% 
% 
% %definição da estrutura que armazenará a aceleração
% f = struct( 'x', imu_measurements.ax, 'y' , imu_measurements.ay , 'z', imu_measurements.az );
% 
% %atualiza a estrutura de distância ângular
% s.x = imu_measurements.wx*T;
% s.y = imu_measurements.wy*T;
% s.z = imu_measurements.wz*T;
% 
% v = sqrt(s.x^2 + s.y^2 + s.z^2);
% 
% 
% %Calcula a matriz de rotação
% A = quaternions2dcm(q);
% 
% %Cria o vetor cuja as medidas sao as variaveis de medidas para a correcao
% %Y = [q_mag_ace(1), q_mag_ace(2), q_mag_ace(3), q_mag_ace(4), 0, 0, 0, 0, 0, A(3,3)*dist_alt, 0]';
% Y = [q_mag_ace(1), q_mag_ace(2), q_mag_ace(3), q_mag_ace(4), 0]';
% 
% %variavel de estado (10x1)
% %X = [q(1), q(2), q(3), q(4), pose_predicted_state.dx_dt, pose_predicted_state.dy_dt, pose_predicted_state.dz_dt, pose_predicted_state.x, pose_predicted_state.y, pose_predicted_state.z]';
% X = [q(1), q(2), q(3), q(4)]';
% 
% 
% 
% %Cria a matriz C em que estão relacionadas as medidas de correção
% C = zeros(5,4);
% 
% %Atitude
% C(1,1) = 1;
% C(2,2) = 1;
% C(3,3) = 1;
% C(4,4) = 1;
% 
% %rz
% %C(10,10) = 1;
% 
% %pseudo observar para que o modulo da matriz de rotação permaneça unitária
% modulo_q = sqrt(q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);
% C(5,1) = - q(1)/modulo_q;
% C(5,2) = - q(2)/modulo_q;
% C(5,3) = - q(3)/modulo_q;
% C(5,4) = - q(4)/modulo_q;
% 
% 
% 
% % Cria a matriz A que está relacionada com o processo. Observar que a
% % matriz A não é a mesmo utilizado pelo Padilha por acreditar que há erros
% % de sinais nas equações por ele escrito.
% A(1,1) = cos(v/2);
% A(2,1) = sin(v/2)*s.x/v;
% A(3,1) = sin(v/2)*s.y/v;
% A(4,1) = sin(v/2)*s.z/v;
% %A(5,1) = (2*q(1)*f.x + 2*f.y - 2*f.z)*T;
% %A(6,1) = (-2*f.x + 2*q(1)*f.y + 2*f.z)*T;
% %A(7,1) = (2*f.x + 2*f.y + 2*q(1)*f.z)*T;
% %A(8,1) = (2*q(1)*f.x + 2*f.y - 2*f.z)*T^2/2;
% %A(9,1) = (-2*f.x + 2*q(1)*f.y + 2*f.z)*T^2/2;
% %A(10,1) = (2*f.x + 2*f.y + 2*q(1)*f.z)*T^2/2;
% 
% 
% A(1,2) = -sin(v/2)*s.x/v;
% A(2,2) = cos(v/2);
% A(3,2) = -sin(v/2)*s.z/v;
% A(4,2) = sin(v/2)*s.y/v;
% %A(5,2) = (2*q(2)*f.x + 2*f.y + 2*f.z)*T;
% %A(6,2) = (2*f.x - 2*q(2)*f.y + 2*f.z)*T;
% %A(7,2) = (2*f.x - 2*f.y - 2*q(2)*f.z)*T;
% %A(8,2) = (2*q(2)*f.x + 2*f.y + 2*f.z)*T^2/2;
% %A(9,2) = (2*f.x - 2*q(2)*f.y + 2*f.z)*T^2/2;
% %A(10,2) = (2*f.x - 2*f.y - 2*q(2)*f.z)*T^2/2;
% 
% A(1,3) = -sin(v/2)*s.y/v;
% A(2,3) = sin(v/2)*s.z/v;
% A(3,3) = cos(v/2);
% A(4,3) = -sin(v/2)*s.x/v;
% %A(5,3) = (-2*q(3)*f.x + 2*f.y - 2*f.z)*T;
% %A(6,3) = (2*f.x + 2*q(3)*f.y + 2*f.z)*T;
% %A(7,3) = (2*f.x - 2*f.y - 2*q(3)*f.z)*T;
% %A(8,3) = (-2*q(3)*f.x + 2*f.y - 2*f.z)*T^2/2;
% %A(9,3) = (2*f.x + 2*q(3)*f.y + 2*f.z)*T^2/2;
% %A(10,3) = (2*f.x - 2*f.y - 2*q(3)*f.z)*T^2/2;
% 
% A(1,4) = -sin(v/2)*s.z/v;
% A(2,4) = -sin(v/2)*s.y/v;
% A(3,4) = sin(v/2)*s.x/v;
% A(4,4) = cos(v/2);
% %A(5,4) = (-2*q(4)*f.x + 2*f.y + 2*f.z)*T;
% %A(6,4) = (-2*f.x - 2*q(4)*f.y + 2*f.z)*T;
% %A(7,4) = (2*f.x + 2*f.y + 2*q(4)*f.z)*T;
% %A(8,4) = (-2*q(4)*f.x + 2*f.y + 2*f.z)*T^2/2;
% %A(9,4) = (-2*f.x - 2*q(4)*f.y + 2*f.z)*T^2/2;
% %A(10,4) = (2*f.x + 2*f.y + 2*q(4)*f.z)*T^2/2;
% 
% %A(1,5) = 0;
% %A(2,5) = 0;
% %A(3,5) = 0;
% %A(4,5) = 0;
% %A(5,5) = 1;
% %A(6,5) = 0;
% %A(7,5) = 0;
% %A(8,5) = T;
% %A(9,5) = 0;
% %A(10,5) = 0;
% 
% %A(1,6) = 0;
% %A(2,6) = 0;
% %A(3,6) = 0;
% %A(4,6) = 0;
% %A(5,6) = 0;
% %A(6,6) = 1;
% %A(7,6) = 0;
% %A(8,6) = 0;
% %A(9,6) = T;
% %A(10,6) = 0;
% 
% %A(1,7) = 0;
% %A(2,7) = 0;
% %A(3,7) = 0;
% %A(4,7) = 0;
% %A(5,7) = 0;
% %A(6,7) = 0;
% %A(7,7) = 1;
% %A(8,7) = 0;
% %A(9,7) = 0;
% %A(10,7) = T;
% 
% %A(1,8) = 0;
% %A(2,8) = 0;
% %A(3,8) = 0;
% %A(4,8) = 0;
% %A(5,8) = 0;
% %A(6,8) = 0;
% %A(7,8) = 0;
% %A(8,8) = 1;
% %A(9,8) = 0;
% %A(10,8) = 0;
% 
% %A(1,9) = 0;
% %A(2,9) = 0;
% %A(3,9) = 0;
% %A(4,9) = 0;
% %A(5,9) = 0;
% %A(6,9) = 0;
% %A(7,9) = 0;
% %A(8,9) = 0;
% %A(9,9) = 1;
% %A(10,9) = 0;
% 
% %A(1,10) = 0;
% %A(2,10) = 0;
% %A(3,10) = 0;
% %A(4,10) = 0;
% %A(5,10) = 0;
% %A(6,10) = 0;
% %A(7,10) = 0;
% %A(8,10) = 0;
% %A(9,10) = 0;
% %A(10,10) = 1;
% 
% 
% 
% %Predição(10x10) da matriz de covariância
% P = A*P*A' + Q;
% 
% 
% %Calculo do ganho de Kalman
% K = P*C'*(C*P*C' + R)';
% 
% 
% %correção do estado
% X = X + K*(Y - C*X);
% 
% %correção da matriz de covariância
% P = (eye(4) - K*C)*P;
% 
% 
% %Monta o novo quaternio de rotação
% q(1) = X(1);
% q(2) = X(2);
% q(3) = X(3);
% q(4) = X(4);
% 
% %Armazena as variaveis de estados corrigidas pelo filtro de Kalman
% [pose_predicted_state.roll,pose_predicted_state.pitch,pose_predicted_state.yaw] = quaternions2euler(q);
% %pose_predicted_state.dx_dt = X(5);
% %pose_predicted_state.dy_dt = X(6);
% %pose_predicted_state.dz_dt = X(7);
% %pose_predicted_state.x = X(8);
% %pose_predicted_state.y = X(9);
% %pose_predicted_state.z = X(10);
