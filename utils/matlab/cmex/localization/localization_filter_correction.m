%Fun��o que a realiza a etapa de corre��o pelo Filtro de Kalman
function  [kf_structure] = localization_filter_correction(filtername, kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)

switch filtername
    case 'ukf2'
        [kf_structure] = localization_filter_correction_ukf2(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    case 'ekf2'
        [kf_structure] = localization_filter_correction_ekf2(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    case 'cekf'
        [kf_structure] = localization_filter_correction_cekf(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    case 'ekf3'
        [kf_structure] = localization_filter_correction_ekf3(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    case 'ekf_decoupled'
        [kf_structure] = localization_filter_correction_ekf_decoupled(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_correction_cekf(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)
%%% corrige usando medidas do algoritmo TRIAD: 
%%% Fun�ao de medi�ao: [eye(4) zeros(4,6)]*X (a medida � dada pelo
%%% estimador TRIAD.
%%% Estimar a medi�ao do quaternion assim como da matriz R usando UT;
if((magnetometermeasure.flagvalidmeasure))
    %%%
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        if kf_structure.flagestimategyrometerbias
            H = [eye(4) zeros(4,12)];
        else
            H = [eye(4) zeros(4,9)];
        end
    else
        if kf_structure.flagestimategyrometerbias
            H = [eye(4) zeros(4,9)];
        else
            H = [eye(4) zeros(4,6)];
        end
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
%%% Fun�ao de medi�ao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
if gpsmeasure.flagvalidpmeasure
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        if kf_structure.flagestimategyrometerbias
            H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3)];
        else
            H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
        end
    else
        if kf_structure.flagestimategyrometerbias
            H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
        else
            H = [zeros(3,4) eye(3,3) zeros(3,3)];
        end
    end
    v = (gpsmeasure.p - H*X);
    K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do GPS:
%%% Fun�ao de medi�ao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
if gpsmeasure.flagvalidvmeasure
    X = kf_structure.X;
    P = kf_structure.P;

    if kf_structure.flagestimateaccelerometerbias
        if kf_structure.flagestimategyrometerbias
            H = [zeros(3,4) zeros(3,3) eyes(3,3) zeros(3,3) zeros(3,3)];
        else
            H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
        end
    else
        if kf_structure.flagestimategyrometerbias
            H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
        else
            H = [zeros(3,4) zeros(3,3) eye(3,3)];
        end
    end
    
    v = (gpsmeasure.v - H*X);
    K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do sonar:
%%% Fun�ao de medi�ao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
if sonarmeasure.flagvalidmeasure
    X = kf_structure.X;
    P = kf_structure.P;
    dHdX1 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(1);
    dHdX2 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(2);
    dHdX3 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(3);
    dHdX4 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(4);
    dHdX7 =  1/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);

    if kf_structure.flagestimateaccelerometerbias
        if kf_structure.flagestimategyrometerbias
            H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0 0 0 0];
        else
            H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0];
        end
    else
        if kf_structure.flagestimategyrometerbias
            H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0];
        else
            H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0];
        end
    end
    
    v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
    K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige a norma do quaternion usando pseudo-medi�oes:
%%% Fun�ao de medi�ao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Fun�ao de medi�ao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
X = kf_structure.X;
P = kf_structure.P;
dHdX1 = 2*X(1);
dHdX2 = 2*X(2);
dHdX3 = 2*X(3);
dHdX4 = 2*X(4);

% TODO: Gyrometer matrices are wrong. Need to fix.
if kf_structure.flagestimateaccelerometerbias
    if kf_structure.flagestimategyrometerbias
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0 0 0 0 0 0 0];
    else
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
    end
else
    if kf_structure.flagestimategyrometerbias
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
    else
        H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
    end
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
function  [kf_structure] = localization_filter_correction_ekf3(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)
%%% corrige usando medidas do algoritmo TRIAD: 
%%% Fun�ao de medi�ao: [eye(4) zeros(4,6)]*X (a medida � dada pelo
%%% estimador TRIAD.
%%% Estimar a medi�ao do quaternion assim como da matriz R usando UT;
if((magnetometermeasure.flagvalidmeasure))
    %%% 
    X = kf_structure.X;
    P = kf_structure.P;
    if kf_structure.flagestimateaccelerometerbias
        H = [eye(4) zeros(4,9)];
    else
        H = [eye(4) zeros(4,6)];
    end
    [Ym,Py] = ConvertedMeasurementTRIAD('ekf3',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
    v = (Ym - H*X);
    K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad_ekf3);
    kf_structure.X = X + K*v;
    X = kf_structure.X;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do GPS:
%%% Fun�ao de medi�ao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
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
%%% Fun�ao de medi�ao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
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
%%% Fun�ao de medi�ao: z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
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

%%% corrige a norma do quaternion usando pseudo-medi�oes:
%%% Fun�ao de medi�ao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Fun�ao de medi�ao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
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
%%% Fun�ao de medi�ao: [eye(4) zeros(4,6)]*X (a medida � dada pelo
%%% estimador TRIAD.
%%% Estimar a medi�ao do quaternion assim como da matriz R usando UT;
if((magnetometermeasure.flagvalidmeasure))
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

%%% corrige usando medidas do GPS:
%%% Fun�ao de medi�ao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
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
%%% Fun�ao de medi�ao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
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
%%% Fun�ao de medi�ao: z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
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

%%% corrige a norma do quaternion usando pseudo-medi�oes:
%%% Fun�ao de medi�ao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Fun�ao de medi�ao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
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
function  [kf_structure] = localization_filter_correction_ekf_decoupled(kf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T)

if magnetometermeasure.flagvalidmeasure
    X = kf_structure.X;
    P = kf_structure.P;

    % H is the measurement function for the magnetometer - in this
    % case, we have:
    % b = A(q)r
    % where b is the actual measurement (data from the magnetometer),
    % A(q) is the DCM written in terms of the quaternions, and r is the
    % reference frame (the local magnetic field in NED coordinates).
    % Therefore: h(X) = A(q)r
    % Here, A(q) = [ q0^2 + q1^2 - q2^2 - q3^2, 2*q1*q2 - 2*q0*q3,                 2*q0*q2 + 2*q1*q3]
    %              [ 2*q0*q3 + 2*q1*q2,         q0^2 - q1^2 + q2^2 - q3^2,         2*q2*q3 - 2*q0*q1]
    %              [ 2*q1*q3 - 2*q0*q2,         2*q0*q1 + 2*q2*q3,                 q0^2 - q1^2 - q2^2 + q3^2]
    % h(X) is the actual magnetometer measurement (normalized) and r is
    % local magnetic field (normalized as well).
    % Reference: model: Unscented Filtering for Spacecraft Attitude
    % Estimation, Crassidis, Markley
    % Rotation matrix: Stardown interial navigation technlogy, page 48
    
    % Calculating the inovation term directly:
    v = [magnetometermeasure.mx; magnetometermeasure.my; magnetometermeasure.mz;]/sqrt(magnetometermeasure.mx.^2+magnetometermeasure.my.^2+magnetometermeasure.mz.^2) - transpose(quaternions2dcm(X(1:4)))*(M/norm(M));
    
    % We also need dh_dX. This is calculated by hand:
    q0 = X(1);
    q1 = X(2);
    q2 = X(3);
    q3 = X(4);
    
    if kf_structure.flagestimateaccelerometerbias
        dh_dX = zeros(3, 13);
    else
        dh_dX = zeros(3, 10);
    end
    
    % First column: derivative in relation to q1
    % 2*m1*q0 + 2*m2*q3 - 2*m3*q2
    % 2*m2*q0 - 2*m1*q3 + 2*m3*q1
    % 2*m1*q2 - 2*m2*q1 + 2*m3*q0
    dh_dX(1, 1) = 2*M(1)/norm(M)*q0 + 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q2;
    dh_dX(2, 1) = 2*M(2)/norm(M)*q0 - 2*M(1)/norm(M)*q3 + 2*M(3)/norm(M)*q1;
    dh_dX(3, 1) = 2*M(1)/norm(M)*q2 - 2*M(2)/norm(M)*q1 + 2*M(3)/norm(M)*q0;
    % Second column: derivative in relation to q2
    % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
    % 2*m1*q2 - 2*m2*q1 + 2*m3*q0
    % 2*m1*q3 - 2*m2*q0 - 2*m3*q1
    dh_dX(1, 2) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
    dh_dX(2, 2) = 2*M(1)/norm(M)*q2 - 2*M(2)/norm(M)*q1 + 2*M(3)/norm(M)*q0;
    dh_dX(3, 2) = 2*M(1)/norm(M)*q3 - 2*M(2)/norm(M)*q0 - 2*M(3)/norm(M)*q1;
    % Third column: derivative in relation to q3
    % 2*m2*q1 - 2*m1*q2 - 2*m3*q0
    % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
    % 2*m1*q0 + 2*m2*q3 - 2*m3*q2
    dh_dX(1, 3) = 2*M(2)/norm(M)*q1 - 2*M(1)/norm(M)*q2 - 2*M(3)/norm(M)*q0;
    dh_dX(2, 3) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
    dh_dX(3, 3) = 2*M(1)/norm(M)*q0 + 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q2;
    % Fourth column: derivative in relation to q4
    % 2*m2*q0 - 2*m1*q3 + 2*m3*q1
    % 2*m3*q2 - 2*m2*q3 - 2*m1*q0
    % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
    dh_dX(1, 4) = 2*M(2)/norm(M)*q0 - 2*M(1)/norm(M)*q3 + 2*M(2)/norm(M)*q1;
    dh_dX(2, 4) = 2*M(3)/norm(M)*q2 - 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q0;
    dh_dX(3, 4) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
    
    % For consistency, call dh_dX = H
    H = dh_dX;
    % Calculate the innovation covariance S
    S = H*P*transpose(H) + kf_structure.ekf_decoupled_magnetometer_R;
    % Kalman gain
    K = P*transpose(H)*inv(S);
    kf_structure.X = X + K*v;
    X = kf_structure.X;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

% Accelerometer Correction - Attitude
% H is the measurement function for the accelerometer - in this
% case, we have the same as the magnetometer. See above for
% details.

X = kf_structure.X;
P = kf_structure.P;

% Calculating the inovation term directly:
G = -G;
v = [imumeasure.ax; imumeasure.ay; imumeasure.az;]/sqrt(imumeasure.ax.^2+imumeasure.ay.^2+imumeasure.az.^2) - transpose(quaternions2dcm(X(1:4)))*(-1.0*G/norm(G));

% We also need dh_dX. This is calculated by hand:
q0 = X(1);
q1 = X(2);
q2 = X(3);
q3 = X(4);

if kf_structure.flagestimateaccelerometerbias
    dh_dX = zeros(3, 13);
else
    dh_dX = zeros(3, 10);
end

% First column: derivative in relation to q1
% 2*m1*q0 + 2*m2*q3 - 2*m3*q2
% 2*m2*q0 - 2*m1*q3 + 2*m3*q1
% 2*m1*q2 - 2*m2*q1 + 2*m3*q0
dh_dX(1, 1) = 2*G(1)/norm(G)*q0 + 2*G(2)/norm(G)*q3 - 2*G(3)/norm(G)*q2;
dh_dX(2, 1) = 2*G(2)/norm(G)*q0 - 2*G(1)/norm(G)*q3 + 2*G(3)/norm(G)*q1;
dh_dX(3, 1) = 2*G(1)/norm(G)*q2 - 2*G(2)/norm(G)*q1 + 2*G(3)/norm(G)*q0;
% Second column: derivative in relation to q2
% 2*m1*q1 + 2*m2*q2 + 2*m3*q3
% 2*m1*q2 - 2*m2*q1 + 2*m3*q0
% 2*m1*q3 - 2*m2*q0 - 2*m3*q1
dh_dX(1, 2) = 2*G(1)/norm(G)*q1 + 2*G(2)/norm(G)*q2 + 2*G(3)/norm(G)*q3;
dh_dX(2, 2) = 2*G(1)/norm(G)*q2 - 2*G(2)/norm(G)*q1 + 2*G(3)/norm(G)*q0;
dh_dX(3, 2) = 2*G(1)/norm(G)*q3 - 2*G(2)/norm(G)*q0 - 2*G(3)/norm(G)*q1;
% Third column: derivative in relation to q3
% 2*m2*q1 - 2*m1*q2 - 2*m3*q0
% 2*m1*q1 + 2*m2*q2 + 2*m3*q3
% 2*m1*q0 + 2*m2*q3 - 2*m3*q2
dh_dX(1, 3) = 2*G(2)/norm(G)*q1 - 2*G(1)/norm(G)*q2 - 2*G(3)/norm(G)*q0;
dh_dX(2, 3) = 2*G(1)/norm(G)*q1 + 2*G(2)/norm(G)*q2 + 2*G(3)/norm(G)*q3;
dh_dX(3, 3) = 2*G(1)/norm(G)*q0 + 2*G(2)/norm(G)*q3 - 2*G(3)/norm(G)*q2;
% Fourth column: derivative in relation to q4
% 2*m2*q0 - 2*m1*q3 + 2*m3*q1
% 2*m3*q2 - 2*m2*q3 - 2*m1*q0
% 2*m1*q1 + 2*m2*q2 + 2*m3*q3
dh_dX(1, 4) = 2*G(2)/norm(G)*q0 - 2*G(1)/norm(G)*q3 + 2*G(2)/norm(G)*q1;
dh_dX(2, 4) = 2*G(3)/norm(G)*q2 - 2*G(2)/norm(G)*q3 - 2*G(3)/norm(G)*q0;
dh_dX(3, 4) = 2*G(1)/norm(G)*q1 + 2*G(2)/norm(G)*q2 + 2*G(3)/norm(G)*q3;

% For consistency, call dh_dX = H
H = dh_dX;
% Calculate the innovation covariance S
accelerometer_R = kf_structure.ekf_decoupled_accelerometer_R*exp(20*abs(norm([imumeasure.ax; imumeasure.ay; imumeasure.az;]) - norm(G)))*eye(3);
if(accelerometer_R(1,1) > 1e3)
    accelerometer_R = 1e3*eye(3);
end

S = H*P*transpose(H) + accelerometer_R;
% Kalman gain
K = P*transpose(H)*inv(S);
kf_structure.X = X + K*v;
X = kf_structure.X;
kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;

%%% corrige usando medidas do GPS:
%%% Fun�ao de medi�ao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
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
%%% Fun�ao de medi�ao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
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
%%% Fun�ao de medi�ao: z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
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

%%% corrige a norma do quaternion usando pseudo-medi�oes:
%%% Fun�ao de medi�ao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Fun�ao de medi�ao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
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
%%% Fun�ao de medi�ao: [eye(4) zeros(4,6)]*X (a medida � dada pelo
%%% estimador TRIAD.
%%% Estimar a medi�ao do quaternion assim como da matriz R usando UT;
if(isstruct(magnetometermeasure) & isstruct(imumeasure))
    %%% 
    X = kf_structure.X;
    P = kf_structure.P;
    H = [eye(4) zeros(4,6)];
    [Ym,Py] = ConvertedMeasurementTRIAD('ukf2',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
    v = (Ym - H*X);
    K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad_ukf);
    kf_structure.X = X + K*v;
    kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
    kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
end

%%% corrige usando medidas do sonar:
%%% Fun�ao de medi�ao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
%%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
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

%%% corrige a norma do quaternion usando pseudo-medi�oes:
%%% Fun�ao de medi�ao: 1 = q0^2 + q1^2 + q2^2 + q3^2
%%% Fun�ao de medi�ao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
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
    case 'ekf3'
        [Ym,Py] = ConvertedMeasurementTRIAD('cekf', X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias);
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
            YSamples(:,i) =  localization_filter_model_gtriad('evaluate', X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
            X = XSamples(:,i+DimX);
            magnetometermeasure.mx = X(1); magnetometermeasure.my = X(2); magnetometermeasure.mz = X(3);
            imumeasure.ax = X(4); imumeasure.ay = X(5); imumeasure.az = X(6); imumeasure.wx = X(7); imumeasure.wy = X(8); imumeasure.wz = X(9);
            YSamples(:,i+DimX) =  localization_filter_model_gtriad('evaluate',X,Px,imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
        end,
        WSamples(1,UtNSamples) = kappa/(DimX+kappa);
        XSamples(:,UtNSamples) = Xm;
        X = XSamples(:,UtNSamples);
        magnetometermeasure.mx = X(1); magnetometermeasure.my = X(2); magnetometermeasure.mz = X(3);
        imumeasure.ax = X(4); imumeasure.ay = X(5); imumeasure.az = X(6); imumeasure.wx = X(7); imumeasure.wy = X(8); imumeasure.wz = X(9);
        YSamples(:,UtNSamples) =  localization_filter_model_gtriad('evaluate',X,Px,imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
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
% %defini��o da estrutura que armazenar� a acelera��o
% f = struct( 'x', imu_measurements.ax, 'y' , imu_measurements.ay , 'z', imu_measurements.az );
% 
% %atualiza a estrutura de dist�ncia �ngular
% s.x = imu_measurements.wx*T;
% s.y = imu_measurements.wy*T;
% s.z = imu_measurements.wz*T;
% 
% v = sqrt(s.x^2 + s.y^2 + s.z^2);
% 
% 
% %Calcula a matriz de rota��o
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
% %Cria a matriz C em que est�o relacionadas as medidas de corre��o
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
% %pseudo observar para que o modulo da matriz de rota��o permane�a unit�ria
% modulo_q = sqrt(q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);
% C(5,1) = - q(1)/modulo_q;
% C(5,2) = - q(2)/modulo_q;
% C(5,3) = - q(3)/modulo_q;
% C(5,4) = - q(4)/modulo_q;
% 
% 
% 
% % Cria a matriz A que est� relacionada com o processo. Observar que a
% % matriz A n�o � a mesmo utilizado pelo Padilha por acreditar que h� erros
% % de sinais nas equa��es por ele escrito.
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
% %Predi��o(10x10) da matriz de covari�ncia
% P = A*P*A' + Q;
% 
% 
% %Calculo do ganho de Kalman
% K = P*C'*(C*P*C' + R)';
% 
% 
% %corre��o do estado
% X = X + K*(Y - C*X);
% 
% %corre��o da matriz de covari�ncia
% P = (eye(4) - K*C)*P;
% 
% 
% %Monta o novo quaternio de rota��o
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