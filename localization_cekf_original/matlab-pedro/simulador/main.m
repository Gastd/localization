%%% Simulador de estimadores de posicionamento tri-dimensional de um
%%% veículo por um Sistema de Navegação Inercial (INS). Os calculos serão
%%% realizados utilizando quaternios.
clear all;
close all; %fecha todas as janelas
clc; %clear screen

%Caminho dos códigos em C
addpath('..\cmex\sensors');
addpath('..\cmex\rotation');
addpath('..\cmex\localization');
addpath('..\cmex\imm');

%escolha da tragetória a ser seguida
trajectory_name = 'trajectory_helix';
%trajectory_name = 'trajectory_hovering';

flagestimatorTRIAD = 0;
flagestimatorRungeKutta = 0;
flagestimatorEKF2 = 0;
flagestimatorCEKF = 0;
flagestimatorHybridEKF =1;

% Flag que determina a possibilidade ou nao de ocorrencia de falhas no 
% magnetometro
flag_MAG_fail = 1;
%Multiplicador do ruido do magnetometro em caso de falha
mag_cov_multiplier = 1000;
var_mag_fault = (3/3)^2;

flagnoise = 1;
flagkf_estimateaccelerometerbias = 1;

%Flag que sinaliza se dados reais devem ser usados (1) ou estes devem ser
%obtidos a partir de uma simulacao (0).
flaguserealdata = 0;
realdatafilename = 'dados/caminhonete_na_unb.mat';

%Se necessario, carrega as leituras de dados reais
if flaguserealdata == 1
    real_data = read_data(realdatafilename);
end

%%% main time parameters:
if flaguserealdata == 1
    T = real_data.Ts;
    t = real_data.t;
    %        t = t(1:200);
else
    T = 0.01; % simulation sampling time
    t = 0:T:20; % time
end

%%% main variables:
%cria as estruturas que conterão as estimativas da posição(e suas derivadas primeiras e segundas ordem)
%%e da atitude do helicoptero, além como o histórico de estados

%cria a estrutura que conterá as posições e as atitudes do helicoptero reais
%dadas pela trajetória

[vehiclestate] = vehicle_getstate(t(1),trajectory_name); % apenas para mostrar a posiçao estimada.
if flaguserealdata == 1
    vehiclestate.x = 0;
    vehiclestate.y = 0;
    vehiclestate.z = 0;
%    [vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw] = quaternions2euler(real_data.q_init2n);
    [rpy] = matquaternions2euler(real_data.q_init2n);
    vehiclestate.roll = rpy(1);
    vehiclestate.pitch = rpy(2);
    vehiclestate.yaw = rpy(3);
else
    %pose_real armazena uma copia da estrutura de dados de pose para cada
    %instante de tempo.
    pose_real = repmat(pose_create,length(t),1);
    
    %Preenche a primeira posicao do vetor de poses.
    pose_real(1) = pose_get_from_vehiclestate(vehiclestate);
end

% Kalman filter pose
flaginitcovariance = 1;

% Matriz de Probabilidades de Transicao (TPM) para a Cadeia de Markov
% que rege a ocorrencia de falhas nas medicoes dos acelerometro e do
% magnetometro.
%
%   m1 = estimacao do TRIAD ocorre normalmente
%   m2 = problemas na estimacao do TRIAD
%
%   TPM(i,j) = P{mk+1=i|mk = j}
TPM = [0.7 0.3; 0.3 0.7];   %Matriz real

%Extrai o numero de modos do sistema a partir da dimensao da TPM
number_modes = size(TPM);
number_modes = number_modes(1);

mode = 1;   %Sensores do TRIAD funcionando
%mode = 2;   %Sensores do TRIAD com problemas

%Estrutura de dados que armazena os diferentes modos em que o sistema opera
if flag_MAG_fail == 1
    imm_modes = zeros(length(t),1);
    imm_modes(1,1) = mode;
end


%%%%%%%%%%%%%%%%%%%% Criacao das estruturas de dados e variaveis de medicao
% vetor gravidade no sistema earth (vetor coluna)
if flaguserealdata == 1
    G = real_data.gn;
else
    G = [0, 0, 9.8182]';
end

%Suposto campo magnetico da terra no sistema inercial (vetor linha). Valor
%máximo em módulo em Bsb 23,83725 micro Tesla.
if flaguserealdata == 1
    M = real_data.mn;
else
    M = [20.80224 -7.8082 -8.63213]';
end

%cria uma estrutura que contera as mediçoes
measurements = struct(  'sonar',repmat(sonar(0, flagnoise),length(t),1),...
    'imu',repmat(imu(0, G, flagnoise, 1),length(t),1),...
    'gps',repmat(gps(0, flagnoise),length(t),1),...
    'magnetometer',repmat(magnetometer(0, M, flagnoise, 1),length(t),1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Variavel auxiliar que tambem contem uma copia da estrutura de
%dados de pose para cada instante de tempo.
aux = repmat(pose_create,length(t),1);

%Dados reais -> preenche com as medicoes dos sensores e seus parametros
if flaguserealdata == 1
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
    [rpy] = matquaternions2euler(real_data.q_init2n);
    aux(1).roll = rpy(1);
    aux(1).pitch = rpy(2);
    aux(1).yaw = rpy(3);    
    aux(1).rollvariance = flaginitcovariance*(4*pi/180)^2;
    aux(1).pitchvariance = flaginitcovariance*(4*pi/180)^2;
    aux(1).yawvariance = flaginitcovariance*(4*pi/180)^2;
    aux(1).xvariance = flaginitcovariance*(1)^2;
    aux(1).yvariance = flaginitcovariance*(1)^2;
    aux(1).zvariance = flaginitcovariance*(1)^2;
    aux(1).dx_dtvariance = flaginitcovariance*(0.1)^2;
    aux(1).dy_dtvariance = flaginitcovariance*(0.1)^2;
    aux(1).dz_dtvariance = flaginitcovariance*(0.1)^2;
else % simulaçao -> utiliza medidas de simulacao
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

if flagestimatorEKF2 == 1
    pose_estimates_ekf2 = aux;
    [ekf2_structure] = localization_filter_init(flagkf_estimateaccelerometerbias);
    [ekf2_structure] = localization_filter_pose2state(pose_estimates_ekf2(1),ekf2_structure,1);
    ekf2_structure.Preset = ekf2_structure.P;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_ekf2 = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_ekf2 = repmat(parameter_estimates_ekf2,length(t),1);
    end
end


if flagestimatorCEKF == 1
    pose_estimates_cekf = aux;
    [cekf_structure] = localization_filter_init(flagkf_estimateaccelerometerbias);
    [cekf_structure] = localization_filter_pose2state(pose_estimates_cekf(1),cekf_structure,1);
    cekf_structure.Preset = cekf_structure.P;
    if flagkf_estimateaccelerometerbias
        parameter_estimates_cekf = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_cekf = repmat(parameter_estimates_cekf,length(t),1);
    end
end

%Estimador de estados hibrido que utiliza dois EKF para compensar
%problemas de leitura no magnetometro.
if flagestimatorHybridEKF == 1
    
%     global H_EKF h_mag
%     
%     syms mx my mz q0 q1 q2 q3 x y z vx vy vz bwx bwy bwz  
    
    %Componentes da funcao nao-linear f = [f1 f2 f3]' que relaciona as
    %medidas do magnetometro ao quaternio de atitude. O vetor constante
    %[mx my mz]' representa o campo magnetico local.
    
%     M_rot = [   q0^2+q1^2-q2^2-q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);
%                 2*(q1*q2+q0*q3) q0^2-q1^2+q2^2-q3^2 2*(q2*q3-q0*q1);
%                 2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) q0^2-q1^2-q2^2+q3^2];
%     
%     f = (M_rot')*([mx my mz]');
    
    
%     f1 = (q0^2+q1^2-q2^2-q3^2)*mx + 2*(q1*q2+q0*q3)*my + 2*(q1*q3-q0*q2)*mz;    
%     f2 = 2*(q1*q2-q0*q3)*mx + (q0^2-q1^2+q2^2-q3^2)*my + 2*(q2*q3+q0*q1)*mz;    
%     f3 = 2*(q1*q3+q0*q2)*mx + 2*(q2*q3-q0*q1)*my + (q0^2-q1^2-q2^2+q3^2)*mz;
    
    %Funcao nao-linear que relaciona a medida do magnetometro ao quaternio
    %de atitude
%     f1 = f(1,1);
%     f2 = f(2,1);
%     f3 = f(3,1);
    
%     h_mag = [f1 f2 f3]';   
%     
%     fprintf('Derivando...');
%     
%     if flagkf_estimateaccelerometerbias == 1
%         %Jacobiana considerando os biases dos girometros
%         H_EKF = [   diff(f1,q0) diff(f1,q1) diff(f1,q2) diff(f1,q3) diff(f1,x) diff(f1,y) diff(f1,z) diff(f1,vx) diff(f1,vy) diff(f1,vz) diff(f1,bwx) diff(f1,bwy) diff(f1,bwz);
%                     diff(f2,q0) diff(f2,q1) diff(f2,q2) diff(f2,q3) diff(f2,x) diff(f2,y) diff(f2,z) diff(f2,vx) diff(f2,vy) diff(f2,vz) diff(f2,bwx) diff(f2,bwy) diff(f2,bwz);
%                     diff(f3,q0) diff(f3,q1) diff(f3,q2) diff(f3,q3) diff(f3,x) diff(f3,y) diff(f3,z) diff(f3,vx) diff(f3,vy) diff(f3,vz) diff(f3,bwx) diff(f3,bwy) diff(f3,bwz);];
%     else
%         %Jacobiana sem os biases dos girometros
%         H_EKF = [   diff(f1,q0) diff(f1,q1) diff(f1,q2) diff(f1,q3) diff(f1,x) diff(f1,y) diff(f1,z) diff(f1,vx) diff(f1,vy) diff(f1,vz);
%                     diff(f2,q0) diff(f2,q1) diff(f2,q2) diff(f2,q3) diff(f2,x) diff(f2,y) diff(f2,z) diff(f2,vx) diff(f2,vy) diff(f2,vz);
%                     diff(f3,q0) diff(f3,q1) diff(f3,q2) diff(f3,q3) diff(f3,x) diff(f3,y) diff(f3,z) diff(f3,vx) diff(f3,vy) diff(f3,vz);];
%     end
%     
%     fprintf('Terminado.\n\n');
       
    % Atribui a pose_estimates_hekf a estrutura de dados que armazena as
    % variaveis de estado para todos os instantes de tempo.
    pose_estimates_hekf = aux;
    
    %Inicializa os parametros dos filtros (numero de variaveis de estado,
    %matrizes de covariancia associadas,etc.)
    [hekf_structure1] = localization_filter_init(flagkf_estimateaccelerometerbias);
    [hekf_structure2] = localization_filter_init(flagkf_estimateaccelerometerbias);
        
    %Traduz a pose do corpo para valores de suas variaveis de estado
    [hekf_structure1] = localization_filter_pose2state(pose_estimates_hekf(1),hekf_structure1,1);
    [hekf_structure2] = localization_filter_pose2state(pose_estimates_hekf(1),hekf_structure2,1);
    
    %A matriz de covariances Preset eh utilizada quando o Filtro de Kalman
    %precisa ser reinicializado (P deixa de ser definida positiva, por
    %exemplo).
    hekf_structure1.Preset = hekf_structure1.P;
    hekf_structure2.Preset = hekf_structure2.P;
    
    %Se necessario, cria a estrutura para armazenar os dados de estimacao
    %de bias dos girometros.
    if flagkf_estimateaccelerometerbias
        parameter_estimates_hekf = struct('bias_ax',0,'bias_ay',0,'bias_az',0);
        parameter_estimates_hekf = repmat(parameter_estimates_hekf,length(t),1);
    end
            
    % Esta estrutura armazena os vetor de estado e matrizes de covariancia
    % gerados como resultado final do IMM.
    hekf_structure_IMM = hekf_structure1;
            
    %Desconhecimento total acerca da matriz inicial
    TPM_est = ones(number_modes,number_modes)*(1/number_modes);           
    %TPM_est = [0.9 0.2;0.1 0.8];
    
    % Matriz dos parametros alpha das distribuicoes de Dirichlet para cada
    % uma das colunas da TPM.
    alpha_matrix = ones(number_modes,number_modes);
    
    %Vetor de probabilidade dos modos
    pModos = ones(number_modes,1)*(1/number_modes);
    
    %Matrizes de saida (C) para os diferentes modos
%     C1 = eye(4);
%     C2 = eye(4);
%     C_matrices = zeros(4,4,2);
%     C_matrices(:,:,1) = C1;
%     C_matrices(:,:,2) = C2;

    %Matrizes de covariancias de medicao (R) para os diferentes modos        
    mag_var_x = measurements.magnetometer(1).mxvariance;
    mag_var_y = measurements.magnetometer(1).myvariance;
    mag_var_z = measurements.magnetometer(1).mzvariance;
    
    %Matriz de covariancias para o magnetometro funcionando normalmente
    R1 = [mag_var_x 0 0; 0 mag_var_y 0; 0 0 mag_var_z];
    
    %Supoe que as medidas podem estar situadas com qualquer valor no
    %intervalo 3sigma   
    R2 = eye(3)*(var_mag_fault);
    
    
%     %Verifica se o magnetometro esta sujeito a falhas
%     if flag_MAG_fail == 1
%         R2 = R1*mag_cov_multiplier;
%     else
%         R2 = R1;
%     end
 
    %Armazena as matrizes de covariancias dos ruidos de medicao do
    %magnetometro.
    %R_matrices = zeros(4,4,2);
    R_matrices(:,:,1) = R1;
    R_matrices(:,:,2) = R2;
    
    %Limpa as variaveis auxiliares
    clear mag_var_x mag_var_y mag_var_z;
    
    %Cria uma estrutura para armazenar simultaneamente os vetores de estado
    %dos dois estimadores.
    dimState = size(hekf_structure1.X);
    state_vectors = zeros(dimState(1),dimState(2),2);
    
    %Cria uma estrutura para armazenar simultaneamente as matrizes de
    %covariancia dos dois estimadores.
    dimCov = size(hekf_structure1.P);
    cov_matrices = zeros(dimCov(1),dimCov(2),2);
    
    %Estrutura de dados para armazenar as probabilidades dos modos
    imm_mode_prob = zeros(number_modes,length(t));
    
    %Estrutura de dados para armazenar as matrizes de covariancia estimadas
    %pela transformada Unscented.
%     Py_1_UT = zeros(4,4,length(t));
%     Py_2_UT = zeros(4,4,length(t));
    
    %Estrutura de dados para armazenar o vetor de verossimilhanca das
    %medidas.
    likelihood_vector_imm = zeros(number_modes,length(t)-1);
    
    %Quaternio do passo anterior de estimacao
%     q_prev = [pose_real(1).q0 pose_real(1).q1 pose_real(1).q2 pose_real(1).q3]';    
end

% 4ª order aproximation pose
if flagestimatorRungeKutta == 1
    pose_estimates_runge_kutta = aux;
end

% TRIAD pose (attitude only)
if flagestimatorTRIAD == 1
    pose_estimates_triad = aux;
end

clear aux;


%%% main view:
%inicia o gráfico do veículo
figure('Name','Grafico do Veiculo','NumberTitle','off'); plot3(0,0,0); hold on; grid on;
%desenha o estado atual do veículo
hvehicle = vehicle_draw(vehiclestate);
if flaguserealdata == 1
    dx = 1000;
    dy = 1000;
    dz = 100;
    axis([vehiclestate.x-dx vehiclestate.x+dx vehiclestate.y-dy vehiclestate.y+dy vehiclestate.z-dz vehiclestate.z+dz]);
else
    %define as partes visíveis dos eixos x, y e z.
    axis([-8 8 -8 8 -1 10]);
end
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
set(gca,'DataAspectRatio',[1 1 1]);

%chama a função que desenhará a seta orientada que caracteriza a attitude
%do veículo
if flaguserealdata == 1
    coordinate_system_draw([vehiclestate.x vehiclestate.y vehiclestate.z]',euler2dcm(vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw),5,'X_e','Y_e','Z_e','m');
else
    coordinate_system_draw([0 0 0]',eye(3),5,'X_e','Y_e','Z_e','m');
end
npreview = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=2:length(t),
    %for n=2:2,
    if rem(n,10)==0
        %disp(sprintf('*** t = %f',t(n)));
        fprintf('*** t = %f\n',t(n));
    end

    
    %Verifica se o magnetometro esta sujeito a falhas
    if flag_MAG_fail == 1     
        
        %Simula as transicoes de uma cadeia de Markov
        %
        %   mode == 1 -> Sensores funcionando corretamente
        %   mode == 2 -> Sensores corrompidos por ruido
        transition_prob = rand(1);            
        for i=1:number_modes
            if transition_prob <= TPM(i,mode)
                mode = i;
                break;
            else
                transition_prob = transition_prob-TPM(i,mode);
            end
        end       
               
%         if(n<(length(t)/3))
%             mode = 1;
%         else            
%             if(n<(length(t)/3)*2)            
%                 mode = 2;
%             else
%                 mode = 1;
%             end
%         end

       

        %Armazena o modo do sistema
        imm_modes(n,1) = mode;
        
        %Multiplicador de incerteza para calculo das matrizes de
        %covariancias.
        %uncertainty_multiplier = 5;
        %uncertainty_multiplier = 1;
        
        %Se o magnetometro esta funcionando corretamente, nao modifica a
        %variancia do ruido.
%         if mode == 1            
%             mag_noise_multiplier = 1;            
%         %Se o magnetometro apresenta falhas, multiplica a variancia do
%         %ruido.
%         else            
%             mag_noise_multiplier = mag_cov_multiplier;
%         end
        
    else
        %Nao modifica a variancia do magnetometro
        mag_noise_multiplier = 1;
    end
    
    % Nao altera a matriz de covariancias dos ruidos das IMU. Este comando
    % esta aqui por motivos de compatibilidade com versoes anteriores do
    % script.
    accel_noise_multiplier = 1;    
    mag_noise_multiplier = 1;
    
    if flaguserealdata
    else
        %%% simulation of the vehicle's motion:
        %adquire a trajetória e armazena
        [vehiclestate] = vehicle_getstate(t(n),trajectory_name);
        pose_real(n) = pose_get_from_vehiclestate(vehiclestate);
        q = quaternions_correctsign([pose_real(n).q0 pose_real(n).q1 pose_real(n).q2 pose_real(n).q3], [pose_real(n-1).q0 pose_real(n-1).q1 pose_real(n-1).q2 pose_real(n-1).q3]);
        pose_real(n).q0 = q(1);    pose_real(n).q1 = q(2);    pose_real(n).q2 = q(3);    pose_real(n).q3 = q(4);
    end

    %%% IMU measurements:
    if flaguserealdata == 1 % variancia default.
        measurements.imu(n).wx = real_data.wx_tilde(n) - real_data.wx_bias;
        measurements.imu(n).wy = real_data.wy_tilde(n) - real_data.wy_bias;
        measurements.imu(n).wz = real_data.wz_tilde(n) - real_data.wz_bias;
        measurements.imu(n).ax = real_data.ax(n);
        measurements.imu(n).ay = real_data.ay(n);
        measurements.imu(n).az = real_data.az(n);
    else       
        [measurements.imu(n)] = imu(vehiclestate, G, flagnoise,accel_noise_multiplier);
    end

    %%% Magnetometer measurements:
    measurements.magnetometer(n).flagvalidmeasure = 1;
    if flaguserealdata == 1 % variancia default.
        measurements.magnetometer(n).mx = real_data.mx(n);
        measurements.magnetometer(n).my = real_data.my(n);
        measurements.magnetometer(n).mz = real_data.mz(n);
        if ((abs(real_data.mx(n))+abs(real_data.my(n))+abs(real_data.mz(n)))==0)
            measurements.magnetometer(n).flagvalidmeasure = 0;
        end
    else
        
        [measurements.magnetometer(n)] = magnetometer(vehiclestate, M, flagnoise,mag_noise_multiplier);
        
        %Gera medicoes simuladas do magnetometro corrompidas por ruido
        if mode == 2                         
            
             %Falha do magnetometro
%             if rand(1)<0.6
%                 measurements.magnetometer(n).mx = 0;
%             else
%                 measurements.magnetometer(n).mx = randn(1)*(var_mag_fault);
%             end
%             
%             if rand(1)<0.6
%                 measurements.magnetometer(n).my = 0;
%             else
%                 measurements.magnetometer(n).my = randn(1)*(var_mag_fault);
%             end
%             
%             if rand(1)<0.6
%                 measurements.magnetometer(n).mz = 0;
%             else
%                 measurements.magnetometer(n).mz = randn(1)*(var_mag_fault);
%             end

            measurements.magnetometer(n).mx = randn(1)*(var_mag_fault);
            measurements.magnetometer(n).my = randn(1)*(var_mag_fault);
            measurements.magnetometer(n).mz = randn(1)*(var_mag_fault);
            
            measurements.magnetometer(n).flagvalidmeasure = 2;
        end
    end

    %%% Sonar measurements:
    if flaguserealdata % variancia default.
    else
        [measurements.sonar(n)] = sonar(vehiclestate, flagnoise);
    end

    %%% GPS measurements:
    measurements.gps(n).flagvalidpmeasure = 1;
    measurements.gps(n).flagvalidvmeasure = 1;
    if flaguserealdata == 1 % variancia default.
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
    if flagestimatorTRIAD == 1
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

    %%% Estimacao baseada no Filtro de Kalman Estendido
    if flagestimatorEKF2 == 1
        % prediçao
%         [ekf2_structure] = mexlocalization('FILTER_EKF2_PREDICTION', ekf2_structure, measurements.imu(n), G, T);
        [ekf2_structure] = localization('FILTER_EKF2_PREDICTION', ekf2_structure, measurements.imu(n), G, T);

        % descomente a proxima linha para correçao por IMU + Magnetometro + GPS + Sonar
        % measurements.sonar(n).flagvalidmeasure = 1;
        % descomente a proxima linha para correçao por IMU + Magnetometro + GPS
        measurements.sonar(n).flagvalidmeasure = 0;
        
        %Altera a flag de medida valida o magnetometro para que a correcao
        %seja feita sem o TRIAD.
        measurements.magnetometer(n).flagvalidmeasure = 2;
        
%         [ekf2_structure] = mexlocalization('FILTER_EKF2_CORRECTION', ekf2_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
        [ekf2_structure] = localization('FILTER_EKF2_CORRECTION', ekf2_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
        
        pose_estimates_ekf2(n) = localization('FILTER_STATE2POSE',pose_estimates_ekf2(n), ekf2_structure,1);
        if flagkf_estimateaccelerometerbias
            parameter_estimates_ekf2(n).bias_ax = ekf2_structure.X(11);
            parameter_estimates_ekf2(n).bias_ay = ekf2_structure.X(12);
            parameter_estimates_ekf2(n).bias_az = ekf2_structure.X(13);
        end
    end

    %%% Estimacao baseada no Filtro de Kalman Estendido Correlato
    if flagestimatorCEKF == 1
        % prediçao
        [cekf_structure] = mexlocalization('FILTER_CEKF_PREDICTION', cekf_structure, measurements.imu(n), G, T);
        
        % descomente a proxima linha para correçao por IMU + Magnetometro + GPS + Sonar
        %[cekf_structure] = localization('FILTER_CEKF_CORRECTION', cekf_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
        % descomente a proxima linha para correçao por IMU + Magnetometro + GPS
        measurements.sonar(n).flagvalidmeasure = 0;
        [cekf_structure] = mexlocalization('FILTER_CEKF_CORRECTION', cekf_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);

        pose_estimates_cekf(n) = localization('FILTER_STATE2POSE',pose_estimates_cekf(n), cekf_structure,1);

        if flagkf_estimateaccelerometerbias
            parameter_estimates_cekf(n).bias_ax = cekf_structure.X(11);
            parameter_estimates_cekf(n).bias_ay = cekf_structure.X(12);
            parameter_estimates_cekf(n).bias_az = cekf_structure.X(13);
        end
    end
    
    
    %%% Estimador hibrido que utiliza dois FKE para tornar-se robusto a
    %%% falhas de medicao do magnetometro
    if flagestimatorHybridEKF == 1 
        
        %   Variaveis de estado do sistema de localizacao
        %
        %   hekf_structure.X(1:4,1) = quaternios de pose (q)
        %   hekf_structure.X(5:7,1) = vetor de posicao linear (p)
        %   hekf_structure.X(8:10,1) = vetor de velocidades lineares (v)
        %   hekf_structure.X(11:13,1) = vetor de biases dos girometros
        %
        %   hekf_structure.P = matriz de covariancias do estado
               
        %Armazena o vetor de probabilidades do último passo
        pModos_ant = pModos;       

        %Updates the mode probability using the Markov TPM
        [pModos_pred] = markov_prob_update(pModos,TPM_est);
        
        %Salva os vetores de estado e as matrizes de covariancia anteriores
        state_vectors(:,:,1) = hekf_structure1.X;
        state_vectors(:,:,2) = hekf_structure2.X;        
        cov_matrices(:,:,1) = hekf_structure1.P;
        cov_matrices(:,:,2) = hekf_structure2.P;
        
        %Mistura de estimativas para geracao das condicoes iniciais para os
        %proximos passos da filtragem
        [mixed_state_vectors,mixed_cov_matrices] = mix_estimates(state_vectors,cov_matrices,pModos,pModos_pred,TPM_est);
        
        %Atribuicao das novas condicoes iniciais
        hekf_structure1.X = mixed_state_vectors(:,:,1);
        hekf_structure2.X = mixed_state_vectors(:,:,2);
        hekf_structure1.P = mixed_cov_matrices(:,:,1);
        hekf_structure2.P = mixed_cov_matrices(:,:,2);       
        
        % Etapa de predicao dos filtros. Os parametros devem ser
        % rigorosamente iguais para os dois filtros, pois ambos
        % compartilham o mesmo modelo de predicao.
%         [hekf_structure1] = mexlocalization('FILTER_EKF2_PREDICTION', hekf_structure1, measurements.imu(n), G, T);
%         [hekf_structure2] = mexlocalization('FILTER_EKF2_PREDICTION', hekf_structure2, measurements.imu(n), G, T);
        [hekf_structure1] = localization('FILTER_EKF2_PREDICTION', hekf_structure1, measurements.imu(n), G, T);
        [hekf_structure2] = localization('FILTER_EKF2_PREDICTION', hekf_structure2, measurements.imu(n), G, T);
             
               
        % Estimativa de posicao dada pelo TRIAD a partir das medicoes do
        % magnetometro e do acelerometro.
%         q_TRIAD = localization_triad(measurements.imu(n), measurements.magnetometer(n), M, G);       
        
        %Corrige o sinal do quaternio estimado pelo TRIAD a partir da
        %predicao de um dos filtros (assume-se que a predicao dos
        %quaternios eh semelhante para todos eles). Isso deve ser feito
        %porque "q" e "-q" sao quaternios que representam a mesma rotacao.       
%         q_TRIAD = quaternions_correctsign(q_TRIAD,q_prev);        
        
        %Atualiza o quaternio anterior
%         q_prev = q_TRIAD;        

        %%% O CALCULO PELA TRANSFORMADA UNSCENTED DA MATRIZ DE COVARIANCIAS
        %%% DO TRIAD NAO EH MAIS NECESSARIO, VISTO QUE AS MEDIDAS DO
        %%% MAGNETOMETRO SAO UTILIZADAS DIRETAMENTE NA CORRECAO DA ATITUDE.
        
%         %Calculo da matriz de covariancias supondo medidas corretas       
%         imu_measure = measurements.imu(n);                
%         mag_measure = measurements.magnetometer(n);        
%         [Ym_1,Py_1] = unscented_transform_triad(imu_measure,mag_measure,M,G,q_TRIAD);        
%         R_matrices(:,:,1) = Py_1;        
%         
%         %Calculo da matriz de covariancias supondo medidas corrompidas        
%         imu_measure = measurements.imu(n);
%         imu_measure.axvariance = imu_measure.axvariance*accel_noise_multiplier*uncertainty_multiplier;
%         imu_measure.ayvariance = imu_measure.ayvariance*accel_noise_multiplier*uncertainty_multiplier;
%         imu_measure.azvariance = imu_measure.azvariance*accel_noise_multiplier*uncertainty_multiplier;
%         
%         mag_measure = measurements.magnetometer(n);
%         mag_measure.mxvariance = mag_measure.mxvariance*mag_noise_multiplier*uncertainty_multiplier;
%         mag_measure.myvariance = mag_measure.myvariance*mag_noise_multiplier*uncertainty_multiplier;
%         mag_measure.mzvariance = mag_measure.mzvariance*mag_noise_multiplier*uncertainty_multiplier;        
%         [Ym_2,Py_2] = unscented_transform_triad(imu_measure,mag_measure,M,G,q_TRIAD);        
%         R_matrices(:,:,2) = Py_2;
%         
%         %Armazena as matrizes para conferencia posterior
%         Py_1_UT(:,:,n-1) = Py_1;
%         Py_2_UT(:,:,n-1) = Py_2;
%         
%         %Testa se a matriz de covariancias Py_1 eh definida positiva
%         if(min(eig(Py_1))<=0)
%             fprintf('Matriz Py_1 mal condicionada\n');
%         end
%         
%         %Testa se a matriz de covariancias Py_2 eh definida positiva
%         if(min(eig(Py_2))<=0)
%             fprintf('Matriz Py_2 mal condicionada\n');
%         end
        
        %Calcula os valores numericos da matriz Jacobiana de medicao H e da
        %funcao de medicao h_est baseada na ultima estimativa de estado.
        [H_1,h_est_1] = subs_EKF(hekf_structure1,M);        
        %[H_2,h_est_2] = subs_EKF(hekf_structure2,M);
        
        %Matriz e funcao de medicao no modo de falha do magnetometro.
        H_2 = zeros(3,13);
        h_est_2 = zeros(3,1);        
                
        %Medidas preditas de cada um dos filtros
        yk_predictions(:,:,1) = h_est_1;
        yk_predictions(:,:,2) = h_est_2;
        
        %Matrizes de covariancias dos termos de inovacao dos EKFs
        P_inov(:,:,1) = H_1*hekf_structure1.P*(H_1') + R_matrices(:,:,1);
        P_inov(:,:,2) = H_2*hekf_structure2.P*(H_2') + R_matrices(:,:,2);          
        
        %Atualiza a probabilidade dos modos a partir das medicoes do 
        %magnetometro e calcula o vetor de verossimilhanca das medidas         
        [pModos,likelihood_vector] = mode_prob_update(pModos_pred,[measurements.magnetometer(n).mx measurements.magnetometer(n).my measurements.magnetometer(n).mz]',yk_predictions,P_inov);
        
        %Armazena o vetor de verossimilhanca para verificao posterior
        likelihood_vector_imm(:,n-1) = likelihood_vector;
        
        % descomente a proxima linha para correçao por IMU + Magnetometro + GPS + Sonar
        %[cekf_structure] = localization('FILTER_CEKF_CORRECTION', cekf_structure, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
        % descomente a proxima linha para correçao por IMU + Magnetometro + GPS
        measurements.sonar(n).flagvalidmeasure = 0;
        
        % Etapa de correcao dos filtros
%         [hekf_structure1] = mexlocalization('FILTER_EKF2_CORRECTION', hekf_structure1, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
%         [hekf_structure2] = mexlocalization('FILTER_EKF2_CORRECTION',
%         hekf_structure2, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);               
        
        %Altera a flag de medida valida o magnetometro para que a correcao
        %seja feita sem o TRIAD.
        measurements.magnetometer(n).flagvalidmeasure = 2;              
        [hekf_structure1] = localization('FILTER_EKF2_CORRECTION', hekf_structure1, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
        
        %Sinaliza que a correcao do segundo filtro deve ser feita
        %assumindo-se que o magnetometro apresenta falhas.
        measurements.magnetometer(n).flagvalidmeasure = 3;        
        [hekf_structure2] = localization('FILTER_EKF2_CORRECTION', hekf_structure2, measurements.gps(n), measurements.imu(n), measurements.magnetometer(n), measurements.sonar(n), M, G, T);
        
        %A saida do IMM eh feita utilizando-se todas as variaveis de estado
        %e suas matrizes de covariancias associadas.
        xk_corrections(:,:,1) = hekf_structure1.X;
        xk_corrections(:,:,2) = hekf_structure2.X;        
        Pk_corrections(:,:,1) = hekf_structure1.P;
        Pk_corrections(:,:,2) = hekf_structure2.P;

        % Vetor de estado e matriz de covariancia gerados pelo algoritmo
        % IMM dados de acordo com as ponderacoes das probabilidades dos
        % modos.      
        [xk_est_IMM,Pk_est_IMM] = imm_output(xk_corrections,Pk_corrections,pModos);
        
        % Armazenamento dos valores de saida no IMM na estrutura de dados
        % padrao do sistema de localizacao.
        hekf_structure_IMM.X = xk_est_IMM;
        hekf_structure_IMM.P = Pk_est_IMM;    
                
        %   Likelihood Dirichlet Estimator
        %
        %   Executa a estimacao da TPM apenas se o vetor de
        %   verossimilhancas contiver entradas nao-nulas. Caso contrario,
        %   problemas de divisao ocorrem.
        if(sum(likelihood_vector)>1e-12)
            [TPM_est,alpha_matrix] = likelihoodDirichletEstimator(likelihood_vector,TPM_est,pModos_ant,alpha_matrix,n-1,[2;2]);
        else
            fprintf('Erro de verossimilhanca = 1: %e 2: %e\n',likelihood_vector(1,1),likelihood_vector(2,1));
        end       
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Este final precisa ser adaptado para que as saidas do filtro
        %%% sejam dadas de acordo com as ponderacoes das probabilidades dos
        %%% modos. -> Acho que ja foi feito
        
        %Traduz as variaveis de pose para as variaveis de estado
        pose_estimates_hekf(n) = localization('FILTER_STATE2POSE',pose_estimates_hekf(n), hekf_structure_IMM,1);
        
        %Se for o caso, armazena as estimativas de bias dos girometros
        if flagkf_estimateaccelerometerbias
            parameter_estimates_hekf(n).bias_ax = hekf_structure_IMM.X(11);
            parameter_estimates_hekf(n).bias_ay = hekf_structure_IMM.X(12);
            parameter_estimates_hekf(n).bias_az = hekf_structure_IMM.X(13);
        end        
        
        %Armazena o vetor de probabilidade dos modos
        imm_mode_prob(:,n) = pModos;       
    end
    

    %%% Pose estimator using 4th order Runge-Kutta:
    if flagestimatorRungeKutta == 1
        pose_estimates_runge_kutta(n) = localization('RUNGEKUTTA',pose_estimates_runge_kutta(n-1), measurements.imu(n), G, T);
    end

    %%% plot update:
    if flaguserealdata == 1
        if real_data.gpsx_n(n)~=0
            vehiclestate.x = real_data.gpsx_n(n);
            vehiclestate.y = real_data.gpsy_n(n);
            vehiclestate.z = real_data.gpsz_n(n);
            hvehicle = vehicle_draw(vehiclestate,hvehicle);
            drawnow;
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
if flaguserealdata == 1
    pose_real = pose_create;
end

% Outros graficos

%Probabilidade dos modos
if flagestimatorHybridEKF == 1
    figure('Name','Probabilidade dos modos','NumberTitle','off');
    plot(t,imm_mode_prob(1,:),'b',t,imm_mode_prob(2,:),'r--');
end




%% graficos
for nestimator=1:4
    switch nestimator
        case 1
            flagplot = flagestimatorEKF2;
            if (flagplot == 0) 
                continue; end
            pose_estimates = pose_estimates_ekf2;
            parameter_estimates = parameter_estimates_ekf2;
            lineformat = 'b-';
            estimatorname = 'EKF2';
        case 2
            flagplot = flagestimatorCEKF;
            if (flagplot == 0) 
                continue; end
            pose_estimates = pose_estimates_cekf;
            parameter_estimates = parameter_estimates_cekf;
            lineformat = 'r-';
            estimatorname = 'CEKF';
        case 3
            flagplot = flagestimatorRungeKutta;
            if (flagplot == 0) 
                continue; end
            pose_estimates = pose_estimates_runge_kutta;
            lineformat = 'm-';
            estimatorname = 'RUNGE';
        case 4
            flagplot = flagestimatorHybridEKF;
            if (flagplot == 0) 
                continue; end
            pose_estimates = pose_estimates_hekf;
            parameter_estimates = parameter_estimates_hekf;
            lineformat = 'g-';
            estimatorname = 'HEKF';
    end

    % grafico da atitude real e suas estimativas em quaternions
    figure('Name','Erro entre a atitude real e suas estimativas em quaternios','NumberTitle','off');

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
    figure('Name','Intervalos 3-sigma das estimativas de atitude em quaternios','NumberTitle','off');

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
    figure('Name','Erro entre a atitude real e suas estimativas em angulos','NumberTitle','off');

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

    % grafico da posição real e suas estimativas
    figure('Name','Erro entre a posição real e suas estimativas','NumberTitle','off');
   
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
    figure('Name','Erro entre a velocidade real e suas estimativas','NumberTitle','off');

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
        figure('Name','Estimativas de bias dos acelerometros','NumberTitle','off');

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

% graficos das mediçoes da IMU
figure('Name','Medicoes da IMU','NumberTitle','off');

subplot(611); plot(t,[measurements.imu.ax],'b'); ylabel('a_x [m/s^2]');
% title('IMU measurements');
subplot(612); plot(t,[measurements.imu.ay],'b'); ylabel('a_y [m/s^2]');
subplot(613); plot(t,[measurements.imu.az],'b'); ylabel('a_z [m/s^2]');
subplot(614); plot(t,[measurements.imu.wx]*180/pi,'b'); ylabel('w_x [deg/s]');
subplot(615); plot(t,[measurements.imu.wy]*180/pi,'b'); ylabel('w_y [deg/s]');
subplot(616); plot(t,[measurements.imu.wz]*180/pi,'b'); ylabel('w_z [deg/s]');
xlabel('t [s]');

% graficos das mediçoes do magnetômetro
figure('Name','Medicoes do magnetometro','NumberTitle','off');

I = find([measurements.magnetometer.flagvalidmeasure]~=0);
data = [measurements.magnetometer.mx];
subplot(311); plot(t(I),data(I),'b'); ylabel('m_x [T]');
% title('Magnetometer measurements');
data = [measurements.magnetometer.my];
subplot(312); plot(t(I),data(I),'b'); ylabel('m_y [T]');
data = [measurements.magnetometer.mz];
subplot(313); plot(t(I),data(I),'b'); ylabel('m_z [T]');
xlabel('t [s]');

% graficos das mediçoes do TRIAD
if flagestimatorTRIAD == 1
    figure('Name','Medicoes do TRIAD','NumberTitle','off');
    
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

%graficos das mediçoes do sonar
figure('Name','Medicoes do sonar','NumberTitle','off');

plot(t,[measurements.sonar.range],'b'); ylabel('r [m]');
%title('Sonar measurements');

%graficos das mediçoes do gps
figure('Name','Medicoes de posicao do GPS','NumberTitle','off');

I = find([measurements.gps.flagvalidpmeasure]~=0);
data = [measurements.gps.p];
subplot(311); plot(t(I),data(1,I),'b.'); ylabel('x [m]');
% title('GPS measurements: position');
subplot(312); plot(t(I),data(2,I),'b.'); ylabel('y [m]');
subplot(313); plot(t(I),data(3,I),'b.'); ylabel('z [m]');
xlabel('t [s]');

%graficos das mediçoes do gps
figure('Name','Medicoes de velocidade do GPS','NumberTitle','off');

I = find([measurements.gps.flagvalidvmeasure]~=0);
data = [measurements.gps.v];
subplot(311); plot(t(I),data(1,I),'b.'); ylabel('v_x [m/s]');
% title('GPS measurements: velocity');
subplot(312); plot(t(I),data(2,I),'b.'); ylabel('v_y [m]');
subplot(313); plot(t(I),data(3,I),'b.'); ylabel('v_z [m]');
xlabel('t [s]');

%grafico 3-D da posiçao
figure('Name','Posicao 3D','NumberTitle','off');

Igps = find([measurements.gps.flagvalidpmeasure]~=0); data = [measurements.gps.p];
plot3(data(1,Igps),data(2,Igps),data(3,Igps),'k.'); hold on;

if flagestimatorEKF2 == 1
    plot3([pose_estimates_ekf2.x],[pose_estimates_ekf2.y],[pose_estimates_ekf2.z],'b');
end

if flagestimatorCEKF == 1
    plot3([pose_estimates_cekf.x],[pose_estimates_cekf.y],[pose_estimates_cekf.z],'r');
end

if flagestimatorHybridEKF == 1
    plot3([pose_estimates_hekf.x],[pose_estimates_hekf.y],[pose_estimates_hekf.z],'g');
end

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
set(gca,'DataAspectRatio',[1 1 1]);

%return;

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
% % graficos das mediçoes da IMU
% figure; %figure 2
% subplot(611); plot(t,[measurements.imu.ax],'b'); ylabel('a_x [m/s^2]');
% title('IMU measurements');
% subplot(612); plot(t,[measurements.imu.ay],'b'); ylabel('a_y [m/s^2]');
% subplot(613); plot(t,[measurements.imu.az],'b'); ylabel('a_z [m/s^2]');
% subplot(614); plot(t,[measurements.imu.wx]*180/pi,'b'); ylabel('w_x [deg/s]');
% subplot(615); plot(t,[measurements.imu.wy]*180/pi,'b'); ylabel('w_y [deg/s]');
% subplot(616); plot(t,[measurements.imu.wz]*180/pi,'b'); ylabel('w_z [deg/s]');
%
% % graficos das mediçoes do magnetômetro e do sonar
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
% % grafico da posição real e suas estimativas
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
