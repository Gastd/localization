function [real_data]  = read_data(filename)

%     Sobre a estrutura dos arquivos, as seguintes variaveis estao disponiveis:
%     - Ts: periodo de amostragem, [s]
%     - N: numero de amostras, ja retiradas as amostras utilizadas no procedimento de inicializacao
%     - t: vetor de tempo, [s]
%     - gravity: magnitude da gravidade local, [m/s^2]
%     - gn: vetor gravidade definido no sistema N ou NED, [m/s^2]
%     - an: vetor forca especifica no sistema N qdo o sistema esta parado (-gn), [m/s^2]
%     - mn: vetor campo magnetico da Terra definido no sistema N [uT]
%     - wx_tilde wy_tilde wz_tilde: medidas DEScalibradas dos girometros, [rad/s]
%     - ax ay az: medidas calibradas dos acelerometros, [m/s^2]
%     - mx my mz: medidas calibradas dos magnetometros, [uT]
%     - gpsx gpsy gpsz: medidas de posicao fornecidas pelo GPS, [m]. Estao no sistema E ou Earth, naturalmente
%     - gpsvx gpsvy gpsvz: medidas de velocidade fornecidas pelo GPS, [m/s]. Tb no sistema E
%     - wx_bias wy_bias wz_bias: estimativas dos biases dos girometros, [rad/s]. Sao obtidas no procedimento de inicializacao, 
%       onde o sistema eh mantido parado
%     - a_init m_init: gravidade e campo magnetico na atitude inicial do sistema. Obtidos na inicializacao
%     - q_init2n: quaternio que transforma coordenadas da atitude inicial do sistema B para o sistema N. Obtido na inicializacao
%     - q_e2n: quaternio que transforma coordenadas do sistema E para o sistema N. Obtido na inicializacao
%     - xn yn zn: coordenadas da origem do sistema N, posicao inicial do sistema, [m]. Obtidas na inicializacao. Lembramos que as 
%       velocidades iniciais sao zero
% 
%     Sobre o procedimento de inicializacao, eu posso mandar os dados referentes a essa parte tb. Oq faço eh, a partir das medidas 
%     iniciais dos acelerometros (a_init) e magnetometros (m_init) do sistema parado, calcular o quaternio do sistema B (ou Body) 
%     inicial para o sistema N, q_b2n (ou q_init2n). Alem disso, acumulo medidas de posicao do GPS para calcular o quaternio 
%     do sistema E para o sistema N, q_e2n, bem como as coordenadas da origem do sistema N representadas no sistema E, (xn,yn,zn).
% 
%     Mando as medidas dos girometros descalibradas ate para lmebrar que para eles nao teve um procedimento de calibracao, 
%     como houve para os acelerometros e magnetometros. Para os girometros, so calculei a media das medidas qdo o sistema estava 
%     parado (na inicializacao) e considerei isso uma estimativa dos biases. Sobre esse assunto, lembro que o modelo dos sensores q adotei foi
% 
%     x_tilde  = x_real*scale + bias + noise 

load(filename);

real_data.N = N;       
real_data.Ts  = Ts;
real_data.a_init = a_init;
real_data.ax = ax;
real_data.ay = ay;
real_data.az = az;
real_data.gravity = gravity;
real_data.m_init = m_init;
real_data.mn = mn;
real_data.mx = mx;
real_data.my = my;
real_data.mz = mz;
real_data.t = t;
real_data.wx_bias = wx_bias;
real_data.wx_tilde = wx_tilde;
real_data.wy_bias = wy_bias;
real_data.wy_tilde = wy_tilde;
real_data.wz_bias = wz_bias;
real_data.wz_tilde = wz_tilde;
real_data.q_init2n = q_init2n;
% real_data.q_e2n = q_e2n;
% real_data.an = an;
real_data.gn = gn;
% real_data.xn = xn;
% real_data.yn = yn;
% real_data.zn = zn;
C_e2n = quaternions2dcm(q_e2n);
p_gps_n = C_e2n*(([gpsx'; gpsy'; gpsz'] - repmat([xn; yn; zn],1,N)));
real_data.gpsx_n = p_gps_n(1,:)';
real_data.gpsy_n = p_gps_n(2,:)';
real_data.gpsz_n = p_gps_n(3,:)';
real_data.gps_validpmeasure = (gpsx ~= 0) & (gpsy ~= 0) & (gpsz ~= 0);
v_gps_n = C_e2n*([gpsvx'; gpsvy'; gpsvz']);
real_data.gpsvx_n = v_gps_n(1,:)';
real_data.gpsvy_n = v_gps_n(2,:)';
real_data.gpsvz_n = v_gps_n(3,:)';
real_data.gps_validvmeasure = (gpsvx ~= 0) & (gpsvy ~= 0) & (gpsvz ~= 0);

