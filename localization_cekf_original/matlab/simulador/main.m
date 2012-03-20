%%% Simulador de estimadores de posicionamento tri-dimensional de um
%%% veï¿½culo por um Sistema de Navegaï¿½ï¿½o Inercial (INS). Os calculos serï¿½o
%%% realizados utilizando quaternios.
clear all;
close all; %fecha todas as janelas
clc; %clear screen

if ispc
    addpath('..\cmex\sensors');
    addpath('..\cmex\rotation');
    addpath('..\cmex\localization');
else
    addpath('../cmex/sensors');
    addpath('../cmex/rotation');
    addpath('../cmex/localization');
end
%escolha da tragetï¿½ria a ser seguida
trajectory_name = 'trajectory_helix';
%trajectory_name = 'trajectory_hovering';

flagestimatorTRIAD = 0;
flagestimatorRungeKutta = 0;
flagestimatorEKF2 = 0;
flagestimatorCEKF = 1;
flagnoise = 1;

flagkf_estimateaccelerometerbias = 1;

flaguserealdata = 0;
if ispc
    realdatafilename = 'dados\caminhonete_na_unb.mat';
else
    realdatafilename = 'dados/caminhonete_na_unb.mat';
end

if flaguserealdata
    real_data = read_data(realdatafilename);
end

%%% main time parameters:
if flaguserealdata
    T = real_data.Ts;
    t = real_data.t;
    %        t = t(1:200);
else
    T = 0.01; % simulation sampling time
    t = [0:T:200]; % time
end

%%% main variables:
%cria as estruturas que conterï¿½o as estimativas da posiï¿½ï¿½o(e suas derivadas primeiras e segundas ordem)
%%e da atitude do helicoptero, alï¿½m como o histï¿½rico de estados

%cria a estrutura que conterï¿½ as posiï¿½ï¿½es e as atitudes do helicoptero reais
%dadas pela trajetï¿½ria
[vehiclestate] = vehicle_getstate(t(1),trajectory_name); % apenas para mostrar a posiï¿½ao estimada.
if flaguserealdata
    vehiclestate.x = 0;
    vehiclestate.y = 0;
    vehiclestate.z = 0;
    [rpy] = quaternions2euler(real_data.q_init2n);
    vehiclestate.roll = rpy(1);
    vehiclestate.pitch = rpy(2);
    vehiclestate.yaw = rpy(3);
else
    pose_real = repmat(pose_create,length(t),1);
    pose_real(1) = pose_get_from_vehiclestate(vehiclestate);
end

% Kalman filter pose
flaginitcovariance = 1;

aux = repmat(pose_create,length(t),1);
if flaguserealdata
    aux(1).x = 0;
    aux(1).y = 0;
    aux(1).z = 0;
    aux(1).dx_dt = 0;
    aux(1).dy_dt = 0;
    aux(1).dz_dt = 0;
    aux(1).q0 = real_data.q_init2n(1);
    aux(1).q1 = real_data.q_init2n(2);
    aux(1).q2 = real_data.q_init2n(3);
    aux(1).q3 = real_data.q_init2n(4);
    [rpy] = quaternions2euler(real_data.q_init2n);
    aux(1).roll = rpy(1);
    aux(1).pitch = rpy(2);
    aux(1).yaw = rpy(3);
    aux(1).rollvariance = flaginitcovariance*(4*pi/180)^2;
    aux(1).pitch = flaginitcovariance*(4*pi/180)^2;
    aux(1).yawvariance = flaginitcovariance*(4*pi/180)^2;
    aux(1).xvariance = flaginitcovariance*(1)^2;
    aux(1).yvariance = flaginitcovariance*(1)^2;
    aux(1).zvariance = flaginitcovariance*(1)^2;
    aux(1).dx_dtvariance = flaginitcovariance*(0.1)^2;
    aux(1).dy_dtvariance = flaginitcovariance*(0.1)^2;
    aux(1).dz_dtvariance = flaginitcovariance*(0.1)^2;
else % simulaï¿½ao
    aux(1) = pose_real(1);
    aux(1).rollvariance = flaginitcovariance*(2*pi/180)^2;
    aux(1).pitchvariance = flaginitcovariance*(2*pi/180)^2;
    aux(1).yawvariance = flaginitcovariance*(2*pi/180)^2;
    aux(1).xvariance = flaginitcovariance*(0.1)^2;
    aux(1).yvariance = flaginitcovariance*(0.1)^2;
    aux(1).zvariance = flaginitcovariance*(0.1)^2;
    aux(1).dx_dtvariance = flaginitcovariance*(0.01)^2;
    aux(1).dy_dtvariance = flaginitcovariance*(0.01)^2;
    aux(1).dz_dtvariance = flaginitcovariance*(0.01)^2;
end
if flagestimatorEKF2
    pose_estimates_ekf2 = aux;
    [ekf2_structure] = localization_filter_init(flagkf_estimateaccelerometerbias);
    [ekf2_structure] = localization_filter_pose2state(pose_estimates_ekf2(1),ekf2_structure,1);
    ekf2_structure.Preset = ekf2_structure.P;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_ekf2 = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_ekf2 = repmat(parameter_estimates_ekf2,length(t),1);
    end
end
if flagestimatorCEKF
    pose_estimates_cekf = aux;
    [cekf_structure] = localization_filter_init(flagkf_estimateaccelerometerbias);
    [cekf_structure] = localization_filter_pose2state(pose_estimates_cekf(1),cekf_structure,1);
    cekf_structure.Preset = cekf_structure.P;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_cekf = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_cekf = repmat(parameter_estimates_cekf,length(t),1);
    end
end

% 4ï¿½ order aproximation pose
if flagestimatorRungeKutta
    pose_estimates_runge_kutta = aux;
end

% TRIAD pose (attitude only)
if flagestimatorTRIAD
    pose_estimates_triad = aux;
end

clear aux;

% vetor gravidade no sistema earth (vetor coluna)
if flaguserealdata
    G = real_data.gn;
else
    G = [0, 0, 9.8182]';
end

%Suposto campo magnetico da terra no sistema inercial (vetor linha). Valor
%mï¿½ximo em mï¿½dulo em Bsb 23,83725 micro Tesla.
if flaguserealdata
    M = real_data.mn;
else
    M = [20.80224 -7.8082 -8.63213]';
end
%cria uma estrutura que conterï¿½ as mediï¿½oes
measurements = struct(  'sonar',repmat(sonar(0, flagnoise),length(t),1),...
    'imu',repmat(imu(0, G, flagnoise),length(t),1),...
    'gps',repmat(gps(0, flagnoise),length(t),1),...
    'magnetometer',repmat(magnetometer(0, M, flagnoise),length(t),1));

%%% main view:
%inicia o grï¿½fico do veï¿½culo
figure(1); plot3(0,0,0); hold on; grid on;
%desenha o estado atual do veï¿½culo
hvehicle = vehicle_draw(vehiclestate);
if flaguserealdata
    dx = 1000;
    dy = 1000;
    dz = 100;
    axis([vehiclestate.x-dx vehiclestate.x+dx vehiclestate.y-dy vehiclestate.y+dy vehiclestate.z-dz vehiclestate.z+dz]);
else
    %define as partes visï¿½veis dos eixos x, y e z.
    axis([-8 8 -8 8 -1 10]);
end
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
set(gca,'DataAspectRatio',[1 1 1]);

%chama a funï¿½ï¿½o que desenharï¿½ a seta orientada que caracteriza a attitude
%do veï¿½culo
if flaguserealdata
    coordinate_system_draw([vehiclestate.x vehiclestate.y vehiclestate.z]',euler2dcm([vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw]),5,'X_e','Y_e','Z_e','m');
else
    coordinate_system_draw([0 0 0]',eye(3),5,'X_e','Y_e','Z_e','m');
end
npreview = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=2:length(t),
    %for n=2:2,
    if rem(n,10)==0
        disp(sprintf('*** t = %f',t(n)));
    end

    if flaguserealdata
    else
        %%% simulation of the vehicle's motion:
        %adquire a trajetï¿½ria e armazena
        [vehiclestate] = vehicle_getstate(t(n),trajectory_name);
        pose_real(n) = pose_get_from_vehiclestate(vehiclestate);
        q = quaternions_correctsign([pose_real(n).q0 pose_real(n).q1 pose_real(n).q2 pose_real(n).q3], [pose_real(n-1).q0 pose_real(n-1).q1 pose_real(n-1).q2 pose_real(n-1).q3]);
        pose_real(n).q0 = q(1);    pose_real(n).q1 = q(2);    pose_real(n).q2 = q(3);    pose_real(n).q3 = q(4);
    end

    %%% IMU measurements:
    if flaguserealdata % variancia default.
        measurements.imu(n).wx = real_data.wx_tilde(n) - real_data.wx_bias;
        measurements.imu(n).wy = real_data.wy_tilde(n) - real_data.wy_bias;
        measurements.imu(n).wz = real_data.wz_tilde(n) - real_data.wz_bias;
        measurements.imu(n).ax = real_data.ax(n);
        measurements.imu(n).ay = real_data.ay(n);
        measurements.imu(n).az = real_data.az(n);
    else
        [measurements.imu(n)] = imu(vehiclestate, G, flagnoise);
    end

    %%% Magnetometer measurements:
    measurements.magnetometer(n).flagvalidmeasure = 1;
    if flaguserealdata % variancia default.
        measurements.magnetometer(n).mx = real_data.mx(n);
        measurements.magnetometer(n).my = real_data.my(n);
        measurements.magnetometer(n).mz = real_data.mz(n);
        if ((abs(real_data.mx(n))+abs(real_data.my(n))+abs(real_data.mz(n)))==0)
            measurements.magnetometer(n).flagvalidmeasure = 0;
        end
    else
        [measurements.magnetometer(n)] = magnetometer(vehiclestate, M, flagnoise);
    end

    %%% Sonar measurements:
    if flaguserealdata % variancia default.
    else
        [measurements.sonar(n)] = sonar(vehiclestate, flagnoise);
    end

    %%% GPS measurements:
    measurements.gps(n).flagvalidpmeasure = 1;
    measurements.gps(n).flagvalidvmeasure = 1;
    if flaguserealdata % variancia default.
        measurements.gps(n).p(1) = real_data.gpsx_n(n);
        measurements.gps(n).p(2) = real_data.gpsy_n(n);
        measurements.gps(n).p(3) = real_data.gpsz_n(n);
        measurements.gps(n).flagvalidpmeasure = real_data.gps_validpmeasure(n);
        measurements.gps(n).v(1) = real_data.gpsvx_n(n);
        measurements.gps(n).v(2) = real_data.gpsvy_n(n);
        measurements.gps(n).v(3) = real_data.gpsvz_n(n);
        measurements.gps(n).flagvalidvmeasure = real_data.gps_validvmeasure(n);
    else
        [measurements.gps(n)] = gps(vehiclestate, flagnoise);
    end

    %%% TRIAD-based attitude estimation
    if flagestimatorTRIAD
        if measurements.magnetometer(n).flagvalidmeasure
            q_previous = [pose_estimates_triad(n-1).q0 pose_estimates_triad(n-1).q1 pose_estimates_triad(n-1).q2 pose_estimates_triad(n-1).q3]';
            q = mexlocalization('TRIAD',measurements.imu(n), measurements.magnetometer(n), M, G, q_previous);
            pose_estimates_triad(n).q0 = q(1);    pose_estimates_triad(n).q1 = q(2);
            pose_estimates_triad(n).q2 = q(3);    pose_estimates_triad(n).q3 = q(4);
            [rpy] = quaternions2euler(q);
            pose_estimates_triad(n).roll = rpy(1);
            pose_estimates_triad(n).pitch = rpy(2);
            pose_estimates_triad(n).yaw = rpy(3);
        end
    end

    %%% Kalman Filter-based estimation
    if flagestimatorEKF2
        % prediï¿½ao
        [ekf2_structure] = mexlocalization('FILTER_EKF2_PREDICTION', ekf2_structure, measurements.imu(n), G, T);

        % descomente a proxima linha para correï¿½ao por IMU + Magnetometro + GPS + Sonar
        % measurements.sonar(n).flagvalidmeasure = 1;
        % descomente a proxima linha para correï¿½ao por IMU + Magnetometro + GPS
        measurements.sonar(n).flagvalidmeasure = 0;
        [ekf2_structure] = mexlocalization('FILTER_EKF2_CORRECTION', ekf2_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);

        pose_estimates_ekf2(n) = localization('FILTER_STATE2POSE',pose_estimates_ekf2(n), ekf2_structure,1);
        if flagkf_estimateaccelerometerbias
            parameter_estimates_ekf2(n).bias_ax = ekf2_structure.X(11);
            parameter_estimates_ekf2(n).bias_ay = ekf2_structure.X(12);
            parameter_estimates_ekf2(n).bias_az = ekf2_structure.X(13);
        end
    end

    if flagestimatorCEKF
        % prediï¿½ao
        [cekf_structure] = mexlocalization('FILTER_CEKF_PREDICTION', cekf_structure, measurements.imu(n), G, T);
        
        % descomente a proxima linha para correï¿½ao por IMU + Magnetometro + GPS + Sonar
        %[cekf_structure] = localization('FILTER_CEKF_CORRECTION', cekf_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
        % descomente a proxima linha para correï¿½ao por IMU + Magnetometro + GPS
        measurements.sonar(n).flagvalidmeasure = 0;
        [cekf_structure] = mexlocalization('FILTER_CEKF_CORRECTION', cekf_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);

        pose_estimates_cekf(n) = localization('FILTER_STATE2POSE',pose_estimates_cekf(n), cekf_structure,1);

        if flagkf_estimateaccelerometerbias
            parameter_estimates_cekf(n).bias_ax = cekf_structure.X(11);
            parameter_estimates_cekf(n).bias_ay = cekf_structure.X(12);
            parameter_estimates_cekf(n).bias_az = cekf_structure.X(13);
        end
    end

    %%% Pose estimator using 4th order Runge-Kutta:
    if flagestimatorRungeKutta
        pose_estimates_runge_kutta(n) = localization('RUNGEKUTTA',pose_estimates_runge_kutta(n-1), measurements.imu(n), G, T);
    end

    %%% plot update:
    if flaguserealdata
        if real_data.gpsx_n(n)~=0
            vehiclestate.x = real_data.gpsx_n(n);
            vehiclestate.y = real_data.gpsy_n(n);
            vehiclestate.z = real_data.gpsz_n(n);
            %hvehicle = vehicle_draw(vehiclestate,hvehicle);
            %drawnow;
        end
    else
        % realiza a cada 200 T.
        if rem(n,ceil(0.2/T))==2
            plot3([pose_real(npreview:n).x],[pose_real(npreview:n).y],[pose_real(npreview:n).z],'r'); %atualiza o grafico
            hvehicle = vehicle_draw(vehiclestate,hvehicle);
            drawnow;
            npreview = n;
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Main Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flagplot3sigma = 0;
sigmacolor = [1 1 1]*0.80;
linewidth = 2;
if flaguserealdata
    pose_real = pose_create;
end

% graficos
for nestimator=1:3
    switch nestimator
        case 1
            flagplot = flagestimatorEKF2;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_ekf2;
            parameter_estimates = parameter_estimates_ekf2;
            lineformat = 'b-';
            estimatorname = 'EKF2';
        case 2
            flagplot = flagestimatorCEKF;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_cekf;
            parameter_estimates = parameter_estimates_cekf;
            lineformat = 'r-';
            estimatorname = 'CEKF';
        case 3
            flagplot = flagestimatorRungeKutta;
            if (flagplot == 0) continue; end
            pose_estimates = pose_estimates_runge_kutta;
            lineformat = 'm-';
            estimatorname = 'CEKF';
    end

    % grafico da atitude real e suas estimativas em quaternions
    figure(2);

    subplot(411); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.q0variance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,([pose_estimates.q0]-[pose_real.q0]),lineformat,'LineWidth',linewidth);
    ylabel('_{q0}');
    %     title(['Error in attitude quaternions: ',estimatorname]);

    subplot(412); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.q1variance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,([pose_estimates.q1]-[pose_real.q1]),lineformat,'LineWidth',linewidth);
    ylabel('{q1}');

    subplot(413); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.q2variance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,([pose_estimates.q2]-[pose_real.q2]),lineformat,'LineWidth',linewidth);
    ylabel('{q2}');

    subplot(414); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.q3variance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,([pose_estimates.q3]-[pose_real.q3]),lineformat,'LineWidth',linewidth);
    ylabel('{q3}');
    xlabel('t [s]');

    % grafico 3-sigma das estimativas de atitude
    figure(3);

    subplot(411); hold on;
    plot(t,3*sqrt([pose_estimates.q0variance]),lineformat,'LineWidth',linewidth);
    ylabel('3\sigma_{q0}');
    %     title(['Estimated 3-\sigma in attitude quaternions: ',estimatorname]);

    subplot(412); hold on;
    plot(t,3*sqrt([pose_estimates.q1variance]),lineformat,'LineWidth',linewidth);
    ylabel('3\sigma_{q1}');

    subplot(413); hold on;
    plot(t,3*sqrt([pose_estimates.q2variance]),lineformat,'LineWidth',linewidth);
    ylabel('3\sigma_{q2}');

    subplot(414); hold on;
    plot(t,3*sqrt([pose_estimates.q3variance]),lineformat,'LineWidth',linewidth);
    ylabel('3\sigma_{q3}');
    xlabel('t [s]');

    % grafico da atitude real e suas estimativas em angulos
    figure(4);

    subplot(311); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.rollvariance])*180/pi;
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    e = [pose_estimates.roll]-[pose_real.roll]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('roll [deg]'); % title(['Error in attitude: ',estimatorname]);

    subplot(312); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.pitchvariance])*180/pi;
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    e = [pose_estimates.pitch]-[pose_real.pitch]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('pitch [deg]');

    subplot(313); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.yawvariance])*180/pi;
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    e = [pose_estimates.yaw]-[pose_real.yaw]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('yaw [deg]'); xlabel('t [s]');

    % grafico da posiï¿½ï¿½o real e suas estimativas
    figure(5);
    subplot(311); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.xvariance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,[pose_estimates.x]-[pose_real.x],lineformat,'LineWidth',linewidth);
    ylabel('x [m]'); % title(['Error in altitude: ',estimatorname]);

    subplot(312); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.yvariance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,[pose_estimates.y]-[pose_real.y],lineformat,'LineWidth',linewidth);
    ylabel('y [m]');

    subplot(313); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.zvariance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,[pose_estimates.z]-[pose_real.z],lineformat,'LineWidth',linewidth);
    ylabel('z [m]'); xlabel('t [s]');

    % grafico da velocidade real e suas estimativas
    figure(6);

    subplot(311); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.dx_dtvariance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,[pose_estimates.dx_dt]-[pose_real.dx_dt],lineformat,'LineWidth',linewidth);
    ylabel('v_x [m/s]'); % title(['Error in velocity: ',estimatorname]);

    subplot(312); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.dy_dtvariance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,[pose_estimates.dy_dt]-[pose_real.dy_dt],lineformat,'LineWidth',linewidth);
    ylabel('v_y [m/s]');

    subplot(313); hold on;
    if(flagplot3sigma)
        sigma = sqrt([pose_estimates.dz_dtvariance]);
        for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
        %plot(t,3*sigma,lineformat,t,-3*sigma,lineformat);
    end
    plot(t,[pose_estimates.dz_dt]-[pose_real.dz_dt],lineformat,'LineWidth',linewidth);
    ylabel('v_z [m/s]'); xlabel('t [s]');

    % grafico das estimativas de bias dos acelerometros
    if flagkf_estimateaccelerometerbias
        figure(7);

        subplot(311); hold on;
        plot(t,[parameter_estimates.bias_ax],lineformat,'LineWidth',linewidth);
        ylabel('b_a_x [m/s^2]'); % title(['Accelerometer bias estimates']);
        subplot(312); hold on;
        plot(t,[parameter_estimates.bias_ay],lineformat,'LineWidth',linewidth);
        ylabel('b_a_y [m/s^2]');
        subplot(313); hold on;
        plot(t,[parameter_estimates.bias_az],lineformat,'LineWidth',linewidth);
        ylabel('b_a_z [m/s^2]');
        xlabel('t [s]');
    end

end

% graficos das mediï¿½oes da IMU
figure;
subplot(611); plot(t,[measurements.imu.ax],'b'); ylabel('a_x [m/s^2]');
% title('IMU measurements');
subplot(612); plot(t,[measurements.imu.ay],'b'); ylabel('a_y [m/s^2]');
subplot(613); plot(t,[measurements.imu.az],'b'); ylabel('a_z [m/s^2]');
subplot(614); plot(t,[measurements.imu.wx]*180/pi,'b'); ylabel('w_x [deg/s]');
subplot(615); plot(t,[measurements.imu.wy]*180/pi,'b'); ylabel('w_y [deg/s]');
subplot(616); plot(t,[measurements.imu.wz]*180/pi,'b'); ylabel('w_z [deg/s]');
xlabel('t [s]');

% graficos das mediï¿½oes do magnetï¿½metro
figure;
I = find([measurements.magnetometer.flagvalidmeasure]~=0);
data = [measurements.magnetometer.mx];
subplot(311); plot(t(I),data(I),'b'); ylabel('m_x [T]');
% title('Magnetometer measurements');
data = [measurements.magnetometer.my];
subplot(312); plot(t(I),data(I),'b'); ylabel('m_y [T]');
data = [measurements.magnetometer.mz];
subplot(313); plot(t(I),data(I),'b'); ylabel('m_z [T]');
xlabel('t [s]');

% graficos das mediï¿½oes do TRIAD
if flagestimatorTRIAD
    figure;
    I = find([measurements.magnetometer.flagvalidmeasure]~=0);
    subplot(311); hold on;
    data = [pose_estimates_triad.roll];
    plot(t(I),data(I)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('roll [deg]'); %title(['Attitude from TRIAD: ',estimatorname]);
    subplot(312); hold on;
    data = [pose_estimates_triad.pitch];
    plot(t(I),data(I)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('pitch [deg]');
    subplot(313); hold on;
    data = [pose_estimates_triad.yaw];
    plot(t(I),data(I)*180/pi,lineformat,'LineWidth',linewidth);
    ylabel('yaw [deg]'); xlabel('t [s]');
end

% graficos das mediï¿½oes do sonar
figure;
plot(t,[measurements.sonar.range],'b'); ylabel('r [m]');
%title('Sonar measurements');

% graficos das mediï¿½oes do gps
figure;
I = find([measurements.gps.flagvalidpmeasure]~=0);
data = [measurements.gps.p];
subplot(311); plot(t(I),data(1,I),'b.'); ylabel('x [m]');
% title('GPS measurements: position');
subplot(312); plot(t(I),data(2,I),'b.'); ylabel('y [m]');
subplot(313); plot(t(I),data(3,I),'b.'); ylabel('z [m]');
xlabel('t [s]');

% graficos das mediï¿½oes do gps
figure;
I = find([measurements.gps.flagvalidvmeasure]~=0);
data = [measurements.gps.v];
subplot(311); plot(t(I),data(1,I),'b.'); ylabel('v_x [m/s]');
% title('GPS measurements: velocity');
subplot(312); plot(t(I),data(2,I),'b.'); ylabel('v_y [m]');
subplot(313); plot(t(I),data(3,I),'b.'); ylabel('v_z [m]');
xlabel('t [s]');

% grafico 3-D da posiï¿½ao
figure;
Igps = find([measurements.gps.flagvalidpmeasure]~=0); data = [measurements.gps.p];
plot3(data(1,Igps),data(2,Igps),data(3,Igps),'k.'); hold on;
if flagestimatorEKF2
    plot3([pose_estimates_ekf2.x],[pose_estimates_ekf2.y],[pose_estimates_ekf2.z],'b');
end
if flagestimatorCEKF
    plot3([pose_estimates_cekf.x],[pose_estimates_cekf.y],[pose_estimates_cekf.z],'r');
end
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
set(gca,'DataAspectRatio',[1 1 1]);

return;
print -depsc -f2 figexp_quaternions.eps;
print -depsc -f3 figexp_quaternions3sigma.eps;
print -depsc -f4 figexp_rollpitchyaw.eps;
print -depsc -f5 figexp_position.eps;
print -depsc -f6 figexp_velocity.eps;
print -depsc -f7 figexp_bias.eps;
print -depsc -f10 figexp_triad.eps;
print -depsc -f12 figexp_positiongps.eps;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Main Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flagplot3sigma = 0;
% sigmacolor = [1 1 1]*0.90;
% formatEKF = 'b-';
% formatCEKF = 'r:';
% linewidth = 2;
%
% % graficos das mediï¿½oes da IMU
% figure; %figure 2
% subplot(611); plot(t,[measurements.imu.ax],'b'); ylabel('a_x [m/s^2]');
% title('IMU measurements');
% subplot(612); plot(t,[measurements.imu.ay],'b'); ylabel('a_y [m/s^2]');
% subplot(613); plot(t,[measurements.imu.az],'b'); ylabel('a_z [m/s^2]');
% subplot(614); plot(t,[measurements.imu.wx]*180/pi,'b'); ylabel('w_x [deg/s]');
% subplot(615); plot(t,[measurements.imu.wy]*180/pi,'b'); ylabel('w_y [deg/s]');
% subplot(616); plot(t,[measurements.imu.wz]*180/pi,'b'); ylabel('w_z [deg/s]');
%
% % graficos das mediï¿½oes do magnetï¿½metro e do sonar
% figure; %figure 3
% subplot(411); plot(t,[measurements.sonar.range],'b'); ylabel('r [m]');
% title('Sonar measurements');
% subplot(412); plot(t,[measurements.magnetometer.mx],'b'); ylabel('m_x [T]');
% title('Magnetometer measurements');
% subplot(413); plot(t,[measurements.magnetometer.my],'b'); ylabel('m_y [T]');
% subplot(414); plot(t,[measurements.magnetometer.mz],'b'); ylabel('m_z [T]');
%
% % grafico da atitude real e suas estimativas em quaternions
% for n
% figure; %figure 4
% subplot(411); hold on;
%     if(flagplot3sigma)
%         if(flagestimatorEKF)
%             sigma = sqrt([pose_estimates_ekf.q0variance]);
%             for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%             plot(t,3*sigma,formatEKF,t,-3*sigma,formatEKF);
%         end
%         if(flagestimatorCEKF)
%             sigma = sqrt([pose_estimates_cekf.q0variance]);
%             for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%             plot(t,3*sigma,formatCEKF,t,-3*sigma,formatCEKF);
%         end
%     end
%     if(flagestimatorRungeKutta) plot(t,([pose_estimates_runge_kutta.q0]-[pose_real.q0]),'b','LineWidth',linewidth); end;
%     if(flagestimatorTRIAD) plot(t,([pose_estimates_triad.q0]-[pose_real.q0]),'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,([pose_estimates_ekf.q0]-[pose_real.q0]),formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,([pose_estimates_cekf.q0]-[pose_real.q0]),formatCEKF,'LineWidth',linewidth); end;
%     ylabel('\epsilon_{q0}');
%     title('Error in attitude quaternions: Runge-Kutta prediction (blue), TRIAD estimation (green), kalman filter (black)');
% subplot(412); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.q1variance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,([pose_estimates_runge_kutta.q1]-[pose_real.q1]),'b','LineWidth',linewidth);  end;
%     if(flagestimatorTRIAD) plot(t,([pose_estimates_triad.q1]-[pose_real.q1]),'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,([pose_estimates_ekf.q1]-[pose_real.q1]),formatEKF,'LineWidth',linewidth);  end;
%     if(flagestimatorCEKF) plot(t,([pose_estimates_cekf.q1]-[pose_real.q1]),formatCEKF,'LineWidth',linewidth);  end;
%     ylabel('\epsilon_{q1}');
% subplot(413); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.q2variance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,([pose_estimates_runge_kutta.q2]-[pose_real.q2]),'b','LineWidth',linewidth);  end;
%     if(flagestimatorTRIAD) plot(t,([pose_estimates_triad.q2]-[pose_real.q2]),'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,([pose_estimates_ekf.q2]-[pose_real.q2]),formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,([pose_estimates_cekf.q2]-[pose_real.q2]),formatCEKF,'LineWidth',linewidth); end;
%     ylabel('\epsilon_{q2}');
% subplot(414); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.q3variance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,([pose_estimates_runge_kutta.q3]-[pose_real.q3]),'b','LineWidth',linewidth); end;
%     if(flagestimatorTRIAD) plot(t,([pose_estimates_triad.q3]-[pose_real.q3]),'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,([pose_estimates_ekf.q3]-[pose_real.q3]),formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,([pose_estimates_cekf.q3]-[pose_real.q3]),formatCEKF,'LineWidth',linewidth); end;
%     ylabel('\epsilon_{q3}');
%     xlabel('t [s]');
%
% figure; subplot(111); hold on;
%     title('Attitude quaternions norm error: Runge-Kutta prediction (blue), TRIAD estimation (green), kalman filter (black)');
%     if(flagestimatorRungeKutta) plot(t,1-sqrt([pose_estimates_runge_kutta.q0].^2+[pose_estimates_runge_kutta.q1].^2+[pose_estimates_runge_kutta.q2].^2+[pose_estimates_runge_kutta.q3].^2),'b','LineWidth',linewidth); end;
%     if(flagestimatorTRIAD) plot(t,1-sqrt([pose_estimates_triad.q0].^2+[pose_estimates_triad.q1].^2+[pose_estimates_triad.q2].^2+[pose_estimates_triad.q3].^2),'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,1-sqrt([pose_estimates_ekf.q0].^2+[pose_estimates_ekf.q1].^2+[pose_estimates_ekf.q2].^2+[pose_estimates_ekf.q3].^2),formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,1-sqrt([pose_estimates_cekf.q0].^2+[pose_estimates_cekf.q1].^2+[pose_estimates_cekf.q2].^2+[pose_estimates_cekf.q3].^2),formatCEKF,'LineWidth',linewidth); end;
%     ylabel('norm');
%     xlabel('t [s]');
%
% % grafico da atitude real e suas estimativas em angulos
% figure; %figure 6
% subplot(311); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.rollvariance])*180/pi;
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     drawnow;
%     if(flagestimatorRungeKutta) e = [pose_estimates_runge_kutta.roll]-[pose_real.roll]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,'b','LineWidth',linewidth); end;
%     if(flagestimatorTRIAD) e = [pose_estimates_triad.roll]-[pose_real.roll]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) e = [pose_estimates_ekf.roll]-[pose_real.roll]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,formatEKF,'LineWidth',linewidth);  end;
%     if(flagestimatorCEKF) e = [pose_estimates_cekf.roll]-[pose_real.roll]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,formatCEKF,'LineWidth',linewidth);  end;
%     ylabel('[degrees]'); title('Error in attitude: Runge-Kutta prediction (blue), TRIAD estimation (green), kalman filter (black)');
% subplot(312); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.pitchvariance])*180/pi;
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) e = [pose_estimates_runge_kutta.pitch]-[pose_real.pitch]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,'b','LineWidth',linewidth); end;
%     if(flagestimatorTRIAD) e = [pose_estimates_triad.pitch]-[pose_real.pitch]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) e = [pose_estimates_ekf.pitch]-[pose_real.pitch]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) e = [pose_estimates_cekf.pitch]-[pose_real.pitch]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,formatCEKF,'LineWidth',linewidth); end;
%     ylabel('[degrees]');
% subplot(313); hold on;
%     if(flagplot3sigma)
%
%         sigma = sqrt([pose_estimates_ekf.yawvariance])*180/pi;
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) e = [pose_estimates_runge_kutta.yaw]-[pose_real.yaw]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,'b','LineWidth',linewidth); end;
%     if(flagestimatorTRIAD) e = [pose_estimates_triad.yaw]-[pose_real.yaw]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,'g','LineWidth',linewidth); end;
%     if(flagestimatorEKF) e = [pose_estimates_ekf.yaw]-[pose_real.yaw]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) e = [pose_estimates_cekf.yaw]-[pose_real.yaw]; e = atan2(sin(e),cos(e)); plot(t,(e)*180/pi,formatCEKF,'LineWidth',linewidth); end;
%     ylabel('[degrees]'); xlabel('t [s]');
%
% % grafico da posiï¿½ï¿½o real e suas estimativas
% figure; % figure 5
% subplot(311); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.xvariance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,[pose_estimates_runge_kutta.x]-[pose_real.x],'b','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,[pose_estimates_ekf.x]-[pose_real.x],formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,[pose_estimates_cekf.x]-[pose_real.x],formatCEKF,'LineWidth',linewidth); end;
%     ylabel('x [m]'); title('Error in altitude: Runge-Kutta prediction (blue), kalman filter (black)')
% subplot(312); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.yvariance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,[pose_estimates_runge_kutta.y]-[pose_real.y],'b','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,[pose_estimates_ekf.y]-[pose_real.y],formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,[pose_estimates_cekf.y]-[pose_real.y],formatCEKF,'LineWidth',linewidth); end;
%     ylabel('y [m]');
% subplot(313); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.zvariance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,[pose_estimates_runge_kutta.z]-[pose_real.z],'b','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,[pose_estimates_ekf.z]-[pose_real.z],formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,[pose_estimates_cekf.z]-[pose_real.z],formatCEKF,'LineWidth',linewidth); end;
%     ylabel('z [m]'); xlabel('t [s]');
%
% % grafico da velocidade real e suas estimativas
% figure; % figure 5
% subplot(311); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.dx_dtvariance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,[pose_estimates_runge_kutta.dx_dt]-[pose_real.dx_dt],'b','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,[pose_estimates_ekf.dx_dt]-[pose_real.dx_dt],formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,[pose_estimates_cekf.dx_dt]-[pose_real.dx_dt],formatCEKF,'LineWidth',linewidth); end;
%     ylabel('v_x [m/s]'); title('Error in velocity: Runge-Kutta prediction (blue), kalman filter (black)')
% subplot(312); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.dy_dtvariance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,[pose_estimates_runge_kutta.dy_dt]-[pose_real.dy_dt],'b','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,[pose_estimates_ekf.dy_dt]-[pose_real.dy_dt],formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,[pose_estimates_cekf.dy_dt]-[pose_real.dy_dt],formatCEKF,'LineWidth',linewidth); end;
%     ylabel('v_y [m/s]');
% subplot(313); hold on;
%     if(flagplot3sigma)
%         sigma = sqrt([pose_estimates_ekf.dz_dtvariance]);
%         for i=2:length(t), h=fill([t(i) t(i-1) t(i-1) t(i)],[3*sigma(i) 3*sigma(i-1) -3*sigma(i-1) -3*sigma(i)],sigmacolor,'EdgeColor',sigmacolor); end
%         plot(t,3*sigma,'m:',t,-3*sigma,'m:');
%     end
%     if(flagestimatorRungeKutta) plot(t,[pose_estimates_runge_kutta.dz_dt]-[pose_real.dz_dt],'b','LineWidth',linewidth); end;
%     if(flagestimatorEKF) plot(t,[pose_estimates_ekf.dz_dt]-[pose_real.dz_dt],formatEKF,'LineWidth',linewidth); end;
%     if(flagestimatorCEKF) plot(t,[pose_estimates_cekf.dz_dt]-[pose_real.dz_dt],formatCEKF,'LineWidth',linewidth); end;
%     ylabel('v_z [m/s]'); xlabel('t [s]');
