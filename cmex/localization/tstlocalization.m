clear all;
close all;
format long;

if ispc
    addpath('..\sensors');
    addpath('..\rotation');
else
    addpath('../sensors');
    addpath('../rotation');
end

flagtestTRIAD = 0;
flagtestFILTER_EKF2_PREDICTION = 0;
flagtestFILTER_EKF2_CORRECTION = 0;
flagtestFILTER_EKF_DECOUPLED_PREDICTION = 0;
flagtestFILTER_EKF_DECOUPLED_CORRECTION = 1;
flagtestFILTER_CEKF_PREDICTION = 0;
flagtestFILTER_CEKF_CORRECTION = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagtestFILTER_CEKF_CORRECTION
    disp('Testando localization(''FILTER_CEKF_CORRECTION'')...');
    
    T = 0.01;
    G = [0, 0, 9.8182]';
    M = [20.80224 -7.8082 -8.63213]';
    flagnoise = 1;
    flagestimateaccelerometerbias = 1;
    
    N = 100;
    for n=1:N
        roll = pi*(2*rand(1)-1);
        pitch = pi/3*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_previous] = mexeuler2quaternions([roll,pitch,yaw]');
        
        roll = pi*(2*rand(1)-1);
        pitch = pi/3*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_real] = mexeuler2quaternions([roll,pitch,yaw]');
        
        vehiclestate.x = 100*randn(1);
        vehiclestate.y = 100*randn(1);
        vehiclestate.z = 100*randn(1);
        vehiclestate.dx_dt = 1*randn(1);
        vehiclestate.dy_dt = 1*randn(1);
        vehiclestate.dz_dt = 1*randn(1);
        vehiclestate.roll = roll;
        vehiclestate.pitch = pitch;
        vehiclestate.yaw = yaw;
        vehiclestate.droll_dt = 0.1*randn(1);
        vehiclestate.dpitch_dt = 0.1*randn(1);
        vehiclestate.dyaw_dt = 0.1*randn(1);
        vehiclestate.d2x_dt2 = 0.1*randn(1);
        vehiclestate.d2y_dt2 = 0.1*randn(1);
        vehiclestate.d2z_dt2 = 0.1*randn(1);
        magnetometermeasure = magnetometer(vehiclestate, M, flagnoise);
        imumeasure = imu(vehiclestate,G,flagnoise);
        gpsmeasure = gps(vehiclestate,flagnoise);
        sonarmeasure = sonar(vehiclestate,flagnoise);
        
        [kf_structure] = localization_filter_init(flagestimateaccelerometerbias);
        kf_structure.X = 0.1*randn(kf_structure.Nstates,1);
        kf_structure.X(1:4) = kf_structure.X(1:4) / norm(kf_structure.X(1:4));
        kf_structure.X(1:4) = kf_structure.X(1:4) + 0.5*randn(4,1);
        A = (randn(kf_structure.Nstates,kf_structure.Nstates));
        kf_structure.P = 0.0001*A*A';
        kf_structure.Preset = eye(kf_structure.Nstates) + A*A;
        
        %        norm(kf_structure.X(1:4))
        
        gpsmeasure.flagvalidpmeasure = 0;
        gpsmeasure.flagvalidvmeasure = 0;
        sonarmeasure.flagvalidmeasure = 0; % +/-
        magnetometermeasure.flagvalidmeasure = 1; %?
        [mexstr,Telapsed(n,1)] = mexlocalization('FILTER_CEKF_CORRECTION',kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T);
        [matstr] = matlocalization('FILTER_CEKF_CORRECTION',{kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T});
        %        [matstr] = matlocalization('FILTER_CEKF_CORRECTION',{kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T});
        
        %return;
        eX = mexstr.X-matstr.X;
        eP = mexstr.P-matstr.P; eP = sum(eP)';
        dX(n) = sqrt(eX'*eX)
        dP(n) = sqrt(eP'*eP)
        
        if (dX(n) > 1e10) | (dP(n) > 1e10)
            disp('Erro para o caso seguinte:');
            roll
            pitch
            yaw
            norm(matstr.X(1:4))
            norm(mexstr.X(1:4))
            disp('Com a seguinte diferen�a:');
            eX
            eP
            return;
        end
    end
    disp('Teste concluido com sucesso!');
    figure; plot(dX);
    title('Erro em X');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(dP);
    title('Erro em P');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(Telapsed*1000);
    title('tempo de execu��o');    xlabel('Numero do teste');    ylabel('Te [ms]');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagtestFILTER_CEKF_PREDICTION
    disp('Testando localization(''FILTER_CEKF_PREDICTION'')...');
    
    T = 0.01;
    G = [0, 0, 9.8182]';
    M = [20.80224 -7.8082 -8.63213]';
    flagnoise = 0;
    flagestimateaccelerometerbias = 0;
    
    N = 1000;
    for n=1:N
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_previous] = mexeuler2quaternions([roll,pitch,yaw]');
        
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_real] = mexeuler2quaternions([roll,pitch,yaw]');
        
        vehiclestate.roll = roll;
        vehiclestate.pitch = pitch;
        vehiclestate.yaw = yaw;
        vehiclestate.droll_dt = 0.1*randn(1);
        vehiclestate.dpitch_dt = 0.1*randn(1);
        vehiclestate.dyaw_dt = 0.1*randn(1);
        vehiclestate.d2x_dt2 = 0.1*randn(1);
        vehiclestate.d2y_dt2 = 0.1*randn(1);
        vehiclestate.d2z_dt2 = 0.1*randn(1);
        magnetometermeasure = magnetometer(vehiclestate, M, flagnoise);
        imumeasure = imu(vehiclestate,G,flagnoise);
        
        [kf_structure] = localization_filter_init(flagestimateaccelerometerbias);
        kf_structure.X = randn(kf_structure.Nstates,1);
        A = (randn(kf_structure.Nstates,kf_structure.Nstates));
        kf_structure.P = A*A';
        kf_structure.Preset = eye(kf_structure.Nstates) + A*A;
        [mexstr,Telapsed(n,1)] = mexlocalization('FILTER_CEKF_PREDICTION',kf_structure,imumeasure,G,T);
        matstr = matlocalization('FILTER_CEKF_PREDICTION',{kf_structure,imumeasure,G,T});
        
        eX = mexstr.X-matstr.X;
        eP = mexstr.P-matstr.P; eP = sum(eP)';
        dX(n) = sqrt(eX'*eX);
        dP(n) = sqrt(eP'*eP);
        if (dX(n) > 1e-10) | (dP(n) > 1e-10)
            disp('Erro para o caso seguinte:');
            roll
            pitch
            yaw
            qmex
            qmat
            disp('Com a seguinte diferen�a:');
            eX
            eP
            return;
        end
    end
    disp('Teste concluido com sucesso!');
    figure; plot(dX);
    title('Erro em X');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(dP);
    title('Erro em P');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(Telapsed*1000);
    title('tempo de execu��o');    xlabel('Numero do teste');    ylabel('Te [ms]');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagtestFILTER_EKF2_CORRECTION
    disp('Testando localization(''FILTER_EKF2_CORRECTION'')...');
    
    T = 0.01;
    G = [0, 0, 9.8182]';
    M = [20.80224 -7.8082 -8.63213]';
    flagnoise = 1;
    flagestimateaccelerometerbias = 1;
    
    N = 1000;
    for n=1:N
        roll = pi*(2*rand(1)-1);
        pitch = pi/3*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_previous] = mexeuler2quaternions([roll,pitch,yaw]');
        
        roll = pi*(2*rand(1)-1);
        pitch = pi/3*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_real] = mexeuler2quaternions([roll,pitch,yaw]');
        
        vehiclestate.x = 100*randn(1);
        vehiclestate.y = 100*randn(1);
        vehiclestate.z = 100*randn(1);
        vehiclestate.dx_dt = 1*randn(1);
        vehiclestate.dy_dt = 1*randn(1);
        vehiclestate.dz_dt = 1*randn(1);
        vehiclestate.roll = roll;
        vehiclestate.pitch = pitch;
        vehiclestate.yaw = yaw;
        vehiclestate.droll_dt = 0.1*randn(1);
        vehiclestate.dpitch_dt = 0.1*randn(1);
        vehiclestate.dyaw_dt = 0.1*randn(1);
        vehiclestate.d2x_dt2 = 0.1*randn(1);
        vehiclestate.d2y_dt2 = 0.1*randn(1);
        vehiclestate.d2z_dt2 = 0.1*randn(1);
        magnetometermeasure = magnetometer(vehiclestate, M, flagnoise);
        imumeasure = imu(vehiclestate,G,flagnoise);
        gpsmeasure = gps(vehiclestate,flagnoise);
        sonarmeasure = sonar(vehiclestate,flagnoise);
        
        [kf_structure] = localization_filter_init(flagestimateaccelerometerbias);
        kf_structure.X = 0.1*randn(kf_structure.Nstates,1);
        kf_structure.X(1:4) = kf_structure.X(1:4) / norm(kf_structure.X(1:4));
        kf_structure.X(1:4) = kf_structure.X(1:4) + 0.5*randn(4,1);
        A = (randn(kf_structure.Nstates,kf_structure.Nstates));
        kf_structure.P = 0.0001*A*A';
        kf_structure.Preset = eye(kf_structure.Nstates) + A*A;
        
        %        norm(kf_structure.X(1:4))
        
        gpsmeasure.flagvalidpmeasure = 0;
        gpsmeasure.flagvalidvmeasure = 0;
        sonarmeasure.flagvalidmeasure = 0; % +/-
        magnetometermeasure.flagvalidmeasure = 1; %?
        [mexstr,Telapsed(n,1)] = mexlocalization('FILTER_EKF2_CORRECTION',kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T);
        [matstr] = matlocalization('FILTER_EKF2_CORRECTION',{kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T});
        %        [matstr] = matlocalization('FILTER_EKF2_CORRECTION',{kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T});
        
        %return;
        eX = mexstr.X-matstr.X;
        eP = mexstr.P-matstr.P; eP = sum(eP)';
        dX(n) = sqrt(eX'*eX);
        dP(n) = sqrt(eP'*eP);
        if (dX(n) > 1e10) | (dP(n) > 1e10)
            disp('Erro para o caso seguinte:');
            roll
            pitch
            yaw
            norm(matstr.X(1:4))
            norm(mexstr.X(1:4))
            disp('Com a seguinte diferen�a:');
            eX
            eP
            return;
        end
    end
    disp('Teste concluido com sucesso!');
    figure; plot(dX);
    title('Erro em X');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(dP);
    title('Erro em P');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(Telapsed*1000);
    title('tempo de execu��o');    xlabel('Numero do teste');    ylabel('Te [ms]');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagtestFILTER_EKF2_PREDICTION
    disp('Testando localization(''FILTER_EKF2_PREDICTION'')...');
    
    T = 0.01;
    G = [0, 0, 9.8182]';
    M = [20.80224 -7.8082 -8.63213]';
    flagnoise = 0;
    flagestimateaccelerometerbias = 0;
    
    N = 1000;
    for n=1:N
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_previous] = mexeuler2quaternions([roll,pitch,yaw]');
        
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_real] = mexeuler2quaternions([roll,pitch,yaw]');
        
        vehiclestate.roll = roll;
        vehiclestate.pitch = pitch;
        vehiclestate.yaw = yaw;
        vehiclestate.droll_dt = 0.1*randn(1);
        vehiclestate.dpitch_dt = 0.1*randn(1);
        vehiclestate.dyaw_dt = 0.1*randn(1);
        vehiclestate.d2x_dt2 = 0.1*randn(1);
        vehiclestate.d2y_dt2 = 0.1*randn(1);
        vehiclestate.d2z_dt2 = 0.1*randn(1);
        magnetometermeasure = magnetometer(vehiclestate, M, flagnoise);
        imumeasure = imu(vehiclestate,G,flagnoise);
        
        [kf_structure] = localization_filter_init(flagestimateaccelerometerbias);
        kf_structure.X = randn(kf_structure.Nstates,1);
        A = (randn(kf_structure.Nstates,kf_structure.Nstates));
        kf_structure.P = A*A';
        kf_structure.Preset = eye(kf_structure.Nstates) + A*A;
        [mexstr,Telapsed(n,1)] = mexlocalization('FILTER_EKF2_PREDICTION',kf_structure,imumeasure,G,T);
        matstr = matlocalization('FILTER_EKF2_PREDICTION',{kf_structure,imumeasure,G,T});
        
        eX = mexstr.X-matstr.X;
        eP = mexstr.P-matstr.P; eP = sum(eP)';
        dX(n) = sqrt(eX'*eX);
        dP(n) = sqrt(eP'*eP);
        if (dX(n) > 1e-10) | (dP(n) > 1e-10)
            disp('Erro para o caso seguinte:');
            roll
            pitch
            yaw
            qmex
            qmat
            disp('Com a seguinte diferen�a:');
            eX
            eP
            return;
        end
    end
    disp('Teste concluido com sucesso!');
    figure; plot(dX);
    title('Erro em X');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(dP);
    title('Erro em P');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(Telapsed*1000);
    title('tempo de execu��o');    xlabel('Numero do teste');    ylabel('Te [ms]');
end

if flagtestFILTER_EKF_DECOUPLED_CORRECTION
    disp('Testando localization(''FILTER_EKF_DECOUPLED_CORRECTION'')...');
    
    T = 0.01;
    G = [0, 0, 9.8182]';
    M = [20.80224 -7.8082 -8.63213]';
    flagnoise = 1;
    flagestimateaccelerometerbias = 1;
    
    N = 1000;
    for n=1:N
        roll = pi*(2*rand(1)-1);
        pitch = pi/3*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_previous] = mexeuler2quaternions([roll,pitch,yaw]');
        
        roll = pi*(2*rand(1)-1);
        pitch = pi/3*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_real] = mexeuler2quaternions([roll,pitch,yaw]');
        
        vehiclestate.x = 100*randn(1);
        vehiclestate.y = 100*randn(1);
        vehiclestate.z = 100*randn(1);
        vehiclestate.dx_dt = 1*randn(1);
        vehiclestate.dy_dt = 1*randn(1);
        vehiclestate.dz_dt = 1*randn(1);
        vehiclestate.roll = roll;
        vehiclestate.pitch = pitch;
        vehiclestate.yaw = yaw;
        vehiclestate.droll_dt = 0.1*randn(1);
        vehiclestate.dpitch_dt = 0.1*randn(1);
        vehiclestate.dyaw_dt = 0.1*randn(1);
        vehiclestate.d2x_dt2 = 0.1*randn(1);
        vehiclestate.d2y_dt2 = 0.1*randn(1);
        vehiclestate.d2z_dt2 = 0.1*randn(1);
        magnetometermeasure = magnetometer(vehiclestate, M, flagnoise);
        imumeasure = imu(vehiclestate,G,flagnoise);
        gpsmeasure = gps(vehiclestate,flagnoise);
        sonarmeasure = sonar(vehiclestate,flagnoise);
        
        [kf_structure] = localization_filter_init(flagestimateaccelerometerbias);
        kf_structure.X = 0.1*randn(kf_structure.Nstates,1);
        kf_structure.X(1:4) = kf_structure.X(1:4) / norm(kf_structure.X(1:4));
        kf_structure.X(1:4) = kf_structure.X(1:4) + 0.5*randn(4,1);
        A = (randn(kf_structure.Nstates,kf_structure.Nstates));
        kf_structure.P = 0.0001*A*A';
        kf_structure.Preset = eye(kf_structure.Nstates) + A*A;
        
        %        norm(kf_structure.X(1:4))
        
        gpsmeasure.flagvalidpmeasure = 0;
        gpsmeasure.flagvalidvmeasure = 0;
        sonarmeasure.flagvalidmeasure = 0; % +/-
        magnetometermeasure.flagvalidmeasure = 0; %?
        [mexstr,Telapsed(n,1)] = mexlocalization('FILTER_EKF_DECOUPLED_CORRECTION',kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T);
        [matstr] = matlocalization('FILTER_EKF_DECOUPLED_CORRECTION',{kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T});
        %        [matstr] = matlocalization('FILTER_EKF2_CORRECTION',{kf_structure,gpsmeasure,imumeasure,magnetometermeasure,sonarmeasure,M,G,T});
        
        %return;
        eX = mexstr.X-matstr.X;
        eP = mexstr.P-matstr.P; eP = sum(eP)';
        dX(n) = sqrt(eX'*eX);
        dP(n) = sqrt(eP'*eP);
        if (dX(n) > 1e10) | (dP(n) > 1e10)
            disp('Erro para o caso seguinte:');
            roll
            pitch
            yaw
            norm(matstr.X(1:4))
            norm(mexstr.X(1:4))
            disp('Com a seguinte diferen�a:');
            eX
            eP
            return;
        end
    end
    disp('Teste concluido com sucesso!');
    figure; plot(dX);
    title('Erro em X');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(dP);
    title('Erro em P');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(Telapsed*1000);
    title('tempo de execu��o');    xlabel('Numero do teste');    ylabel('Te [ms]');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagtestFILTER_EKF_DECOUPLED_PREDICTION
    disp('Testando localization(''FILTER_EKF_DECOUPLED_PREDICTION'')...');
    
    T = 0.01;
    G = [0, 0, 9.8182]';
    M = [20.80224 -7.8082 -8.63213]';
    flagnoise = 0;
    flagestimateaccelerometerbias = 0;
    
    N = 1000;
    for n=1:N
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_previous] = mexeuler2quaternions([roll,pitch,yaw]');
        
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_real] = mexeuler2quaternions([roll,pitch,yaw]');
        
        vehiclestate.roll = roll;
        vehiclestate.pitch = pitch;
        vehiclestate.yaw = yaw;
        vehiclestate.droll_dt = 0.1*randn(1);
        vehiclestate.dpitch_dt = 0.1*randn(1);
        vehiclestate.dyaw_dt = 0.1*randn(1);
        vehiclestate.d2x_dt2 = 0.1*randn(1);
        vehiclestate.d2y_dt2 = 0.1*randn(1);
        vehiclestate.d2z_dt2 = 0.1*randn(1);
        magnetometermeasure = magnetometer(vehiclestate, M, flagnoise);
        imumeasure = imu(vehiclestate,G,flagnoise);
        
        [kf_structure] = localization_filter_init(flagestimateaccelerometerbias);
        kf_structure.X = randn(kf_structure.Nstates,1);
        A = (randn(kf_structure.Nstates,kf_structure.Nstates));
        kf_structure.P = A*A';
        kf_structure.Preset = eye(kf_structure.Nstates) + A*A;
        [mexstr,Telapsed(n,1)] = mexlocalization('FILTER_EKF_DECOUPLED_PREDICTION',kf_structure,imumeasure,G,T);
        matstr = matlocalization('FILTER_EKF_DECOUPLED_PREDICTION',{kf_structure,imumeasure,G,T});
        
        eX = mexstr.X-matstr.X;
        eP = mexstr.P-matstr.P; eP = sum(eP)';
        dX(n) = sqrt(eX'*eX);
        dP(n) = sqrt(eP'*eP);
        if (dX(n) > 1e-10) | (dP(n) > 1e-10)
            disp('Erro para o caso seguinte:');
            roll
            pitch
            yaw
            disp('Com a seguinte diferen�a:');
            eX
            eP
            return;
        end
    end
    disp('Teste concluido com sucesso!');
    figure; plot(dX);
    title('Erro em X');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(dP);
    title('Erro em P');    xlabel('Numero do teste');    ylabel('Norma do erro');
    figure; plot(Telapsed*1000);
    title('tempo de execu��o');    xlabel('Numero do teste');    ylabel('Te [ms]');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flagtestTRIAD
    disp('Testando localization(''TRIAD'')...');
    
    G = [0, 0, 9.8182]';
    M = [20.80224 -7.8082 -8.63213]';
    flagnoise = 0;
    
    N = 1000;
    for n=1:N
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_previous] = mexeuler2quaternions([roll,pitch,yaw]');
        
        roll = pi*(2*rand(1)-1);
        pitch = pi/2*(2*rand(1)-1);
        yaw = pi*(2*rand(1)-1);
        [q_real] = mexeuler2quaternions([roll,pitch,yaw]');
        
        vehiclestate.roll = roll;
        vehiclestate.pitch = pitch;
        vehiclestate.yaw = yaw;
        vehiclestate.droll_dt = 0.1*randn(1);
        vehiclestate.dpitch_dt = 0.1*randn(1);
        vehiclestate.dyaw_dt = 0.1*randn(1);
        vehiclestate.d2x_dt2 = 0.1*randn(1);
        vehiclestate.d2y_dt2 = 0.1*randn(1);
        vehiclestate.d2z_dt2 = 0.1*randn(1);
        magnetometermeasure = magnetometer(vehiclestate, M, flagnoise);
        imumeasure = imu(vehiclestate,G,flagnoise);
        
        [qmex,Telapsed(n,1)]  = mexlocalization('TRIAD',imumeasure,magnetometermeasure,M,G,q_previous);
        qmat = matlocalization('TRIAD',{imumeasure,magnetometermeasure,M,G,q_previous});
        
        e = qmex-qmat;
        d(n) = sqrt(e'*e);
        if d(n) > 1e-10
            disp('Erro para o caso seguinte:');
            roll
            pitch
            yaw
            qmex
            qmat
            disp('Com a seguinte diferen�a:');
            e
            return;
        end
    end
    disp('Teste concluido com sucesso!');
    plot(d);
    xlabel('Numero do teste');
    ylabel('Norma do erro');
    figure; plot(Telapsed*1000);
    title('tempo de execu��o');    xlabel('Numero do teste');    ylabel('Te [ms]');
end
