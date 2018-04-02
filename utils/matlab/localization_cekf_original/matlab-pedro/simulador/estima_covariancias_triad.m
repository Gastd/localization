%   Este script estima a matriz de covariancia dos quaternios gerados pelo
%   metodo TRIAD por meio de simulacoes de Monte Carlo.

clear all;
clc;

%Caminho dos códigos em C
addpath('..\cmex\sensors');
addpath('..\cmex\rotation');
addpath('..\cmex\localization');
addpath('..\cmex\imm');

%escolha da tragetória a ser seguida
trajectory_name = 'trajectory_helix';
%trajectory_name = 'trajectory_hovering';

% Covariance flag
flaginitcovariance = 1;
%Noise flag
flagnoise = 1;

T = 0.01; % simulation sampling time
t = 0:T:10; % time

%Numero de amostras para estimar a matriz de covariancias do TRIAD
num_samples = 5000;

%Multiplicadores das covariancias experimentais dos sensores
accel_noise_multiplier = 1;
mag_noise_multiplier = 1;

%Obtem o valor das variaveis de estado do corpo para uma dada trajetoria e
%em um instante de tempo definido.
[vehiclestate] = vehicle_getstate(t(1),trajectory_name); % apenas para mostrar a posiçao estimada.

%pose_real armazena uma copia da estrutura de dados de pose para cada
%instante de tempo.
pose_real = repmat(pose_create,length(t),1);
%Preenche a primeira posicao do vetor de poses.
pose_real(1) = pose_get_from_vehiclestate(vehiclestate);

%Variavel auxiliar que tambem contem uma copia da estrutura de
%dados de pose para cada instante de tempo.
aux = repmat(pose_create,length(t),1);
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

%Estrutura que armazena as variaveis de pose estimadas pelo TRIAD
pose_estimates_triad = aux;

%Vetor de campo magnetico em Brasilia
M = [20.80224 -7.8082 -8.63213]';

%Vetor de campo gravitacional local
G = [0, 0, 9.8182]';

%cria uma estrutura que conterá as mediçoes
measurements = struct(  'sonar',repmat(sonar(0, flagnoise),length(t),1),...
    'imu',repmat(imu(0, G, flagnoise, 1),length(t),1),...
    'gps',repmat(gps(0, flagnoise),length(t),1),...
    'magnetometer',repmat(magnetometer(0, M, flagnoise, 1),length(t),1));

%Gera as amostras do TRIAD
for n=2:num_samples
           
    fprintf('*** n = %d\n',n);
   
    % Gera uma nova medida da IMU
    [measurements.imu(n)] = imu(vehiclestate, G, flagnoise,accel_noise_multiplier);
    
    % Gera uma nova medida do magnetometro
    measurements.magnetometer(n).flagvalidmeasure = 1;
    [measurements.magnetometer(n)] = magnetometer(vehiclestate, M, flagnoise,mag_noise_multiplier);
    
    % Estima a pose do veiculo por meio do algoritmo TRIAD
    q_previous = [pose_estimates_triad(n-1).q0 pose_estimates_triad(n-1).q1 pose_estimates_triad(n-1).q2 pose_estimates_triad(n-1).q3]';
    q = mexlocalization('TRIAD',measurements.imu(n), measurements.magnetometer(n), M, G, q_previous);
 
    %Armazena a estimativa
    pose_estimates_triad(n).q0 = q(1);    pose_estimates_triad(n).q1 = q(2);
    pose_estimates_triad(n).q2 = q(3);    pose_estimates_triad(n).q3 = q(4);    
end
 
%Calcula o quaternio medio
q0_mean = 0;
q1_mean = 0;
q2_mean = 0;
q3_mean = 0;

for n=1:num_samples
    q0_mean = q0_mean+ pose_estimates_triad(n).q0;
    q1_mean = q1_mean+ pose_estimates_triad(n).q1;
    q2_mean = q2_mean+ pose_estimates_triad(n).q2;
    q3_mean = q3_mean+ pose_estimates_triad(n).q3; 
end

%Quaternio medio
q_mean = ([q0_mean q1_mean q2_mean q3_mean]')*(1/n);

%Calcula a matriz de covariancias
R = zeros(4,4);

for n=1:num_samples
    
    q_i = [pose_estimates_triad(n).q0 pose_estimates_triad(n).q1 pose_estimates_triad(n).q2 pose_estimates_triad(n).q3]';    
    R = R + (q_i-q_mean)*((q_i-q_mean)');
end

%Matriz de covariancias estimada
R = R*(1/(n-1));
 
 
 
 
 
 
 
 
 
 
 
 