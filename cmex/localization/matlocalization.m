function result = matlocalization(proceedurestrg, varargincell)

% Algoritmo TRIAD Improved
if strcmp(proceedurestrg,'TRIAD')
    imumeasure = varargincell{1};
    magnetometermeasure = varargincell{2};
    M = varargincell{3};
    G = varargincell{4};
    q_previous = varargincell{5};

    q = localization_triad(imumeasure, magnetometermeasure, M, G, q_previous);
    result = q;
end

% Algoritmo FKE aplicado � arquitetura correlata: etapa de predi��o
if strcmp(proceedurestrg,'FILTER_EKF2_PREDICTION')
    ekf_structure = varargincell{1};
    imumeasure = varargincell{2};
    G = varargincell{3};
    T = varargincell{4};

    [ekf_structure] = localization_filter_prediction('ekf2', ekf_structure, imumeasure, G, T);
    result = ekf_structure;
end

% Algoritmo FKEC aplicado � arquitetura correlata: etapa de predi��o
if strcmp(proceedurestrg,'FILTER_CEKF_PREDICTION')
    cekf_structure = varargincell{1};
    imumeasure = varargincell{2};
    G = varargincell{3};
    T = varargincell{4};

    [cekf_structure] = localization_filter_prediction('cekf', cekf_structure, imumeasure, G, T);
    result = cekf_structure;
end

if strcmp(proceedurestrg,'FILTER_UKF2_PREDICTION')
    ukf2_structure = varargincell{1};
    imumeasure = varargincell{2};
    G = varargincell{3};
    T = varargincell{4};

    [ukf2_structure] = localization_filter_prediction('ukf2', ukf2_structure, imumeasure, G, T);
    result = ukf2_structure;
end

% Algoritmo FKE aplicado � arquitetura correlata: etapa de corre��o
if strcmp(proceedurestrg,'FILTER_EKF2_CORRECTION')
    ekf_structure = varargincell{1};
    gpsmeasure = varargincell{2};
    imumeasure = varargincell{3};
    magnetometermeasure = varargincell{4};
    sonarmeasure = varargincell{5};
    M = varargincell{6};
    G = varargincell{7};
    T = varargincell{8};

    [ekf_structure] = localization_filter_correction('ekf2', ekf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    result = ekf_structure;
end

% Algoritmo FKEC aplicado � arquitetura correlata: etapa de corre��o
if strcmp(proceedurestrg,'FILTER_CEKF_CORRECTION')
    cekf_structure = varargincell{1};
    gpsmeasure = varargincell{2};
    imumeasure = varargincell{3};
    magnetometermeasure = varargincell{4};
    sonarmeasure = varargincell{5};
    M = varargincell{6};
    G = varargincell{7};
    T = varargincell{8};

    [cekf_structure] = localization_filter_correction('cekf', cekf_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    result = cekf_structure;
end

if strcmp(proceedurestrg,'FILTER_UKF2_CORRECTION')
    ukf2_structure = varargincell{1};
    gpsmeasure = varargincell{2};
    imumeasure = varargincell{3};
    magnetometermeasure = varargincell{4};
    sonarmeasure = varargincell{5};
    M = varargincell{6};
    G = varargincell{7};
    T = varargincell{8};

    [ukf2_structure] = localization_filter_correction('ukf2', ukf2_structure, gpsmeasure, imumeasure, magnetometermeasure, sonarmeasure, M, G, T);
    result = ukf2_structure;
end

if strcmp(proceedurestrg,'FILTER_STATE2POSE')
    pose = varargincell{1};
    ekf_structure = varargincell{2};
    flagpropagateuncertainty = varargincell{3};

    pose = localization_filter_state2pose(pose, ekf_structure,flagpropagateuncertainty);
    result = pose;
end

% Algoritmo RUNGEKUTTA para integra��o do modelo cinem�tico
if strcmp(proceedurestrg,'RUNGEKUTTA')
    pose_estimates_runge_kutta_previous = varargincell{1};
    imumeasure = varargincell{2};
    G = varargincell{3};
    T = varargincell{4};

    pose_estimates_runge_kutta_current = localization_rungekutta(pose_estimates_runge_kutta_previous, imumeasure, G, T);
    result = pose_estimates_runge_kutta_current;
end
% 
% 
