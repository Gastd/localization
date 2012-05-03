%Fun��o que faz a predi��o do filtro de kalman atravez de aproxima��es de
%primeira ordem
function  [kf_structure] = localization_filter_prediction(filtername, kf_structure, imumeasure, G, T)

switch filtername
    case 'ekf2'
        [kf_structure] = localization_filter_prediction_ekf2(kf_structure, imumeasure, G, T);
    case 'cekf'
        [kf_structure] = localization_filter_prediction_cekf(kf_structure, imumeasure, G, T);
    case 'ekf3'
        [kf_structure] = localization_filter_prediction_ekf3(kf_structure, imumeasure, G, T);
    case 'ukf2'
        [kf_structure] = localization_filter_prediction_ukf2(kf_structure, imumeasure, G, T);
    case 'ekf_decoupled'
        [kf_structure] = localization_filter_prediction_ekf_decoupled(kf_structure, imumeasure, G, T);
end

if(sum(eig(kf_structure.P)<0)>0)
    % filter reset:
    disp('*** Warning: Kalman filter reset');
    kf_structure.P = kf_structure.Preset;
end

% if(sum(eig(kf_structure.P)<0)>0)
%     A = kf_structure.Q;
%     B = Py;
%     save matrizes A B;
%     disp('matrices saved')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_prediction_cekf(kf_structure, imumeasure, G, T)
%Etapa de predi��o usando a unscented transform x(k + 1) = f (x(k), u(k + 1))
%     kf_structure.X = f(kf_structure.X,...
%         [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz],...
%         G,T);
%     return;

x_previous = kf_structure.X;
P_previous = kf_structure.P;
u_imu = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
P_u_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
K_previous = kf_structure.S_previous_cekf * inv(kf_structure.R_previous_cekf);
kf_structure.X = localization_filter_model_f('evaluate',x_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias) + K_previous*kf_structure.innovation_previous_cekf;
df_dx = localization_filter_model_f('df_dx',x_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
df_du_imu = localization_filter_model_f('df_du',x_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
Qtilde = kf_structure.Q_cekf + df_du_imu*P_u_imu*df_du_imu'; 
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
kf_structure.P = (df_dx-K_previous*H)*P_previous*(df_dx-K_previous*H)' + Qtilde - K_previous*kf_structure.R_previous_cekf*K_previous';
kf_structure.A_imu_P_imu_cekf = df_du_imu*P_u_imu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_prediction_ekf3(kf_structure, imumeasure, G, T)
%Etapa de predi��o usando a unscented transform x(k + 1) = f (x(k), u(k + 1))
%     kf_structure.X = f(kf_structure.X,...
%         [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz],...
%         G,T);
%     return;

x_previous = kf_structure.X;
P_previous = kf_structure.P;
u = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
Pu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
kf_structure.X = localization_filter_model_f('evaluate',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
df_dx = localization_filter_model_f('df_dx',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
df_du = localization_filter_model_f('df_du',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
kf_structure.P = df_dx*P_previous*df_dx' + df_du*Pu*df_du' + kf_structure.Q_ekf3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_prediction_ekf2(kf_structure, imumeasure, G, T)
%Etapa de predi��o usando a unscented transform x(k + 1) = f (x(k), u(k + 1))
%     kf_structure.X = f(kf_structure.X,...
%         [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz],...
%         G,T);
%     return;

x_previous = kf_structure.X;
P_previous = kf_structure.P;
u = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
Pu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
kf_structure.X = localization_filter_model_f('evaluate',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
df_dx = localization_filter_model_f('df_dx',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
df_du = localization_filter_model_f('df_du',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
kf_structure.P = df_dx*P_previous*df_dx' + df_du*Pu*df_du' + kf_structure.Q_ekf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_prediction_ekf_decoupled(kf_structure, imumeasure, G, T)
%Etapa de predi��o usando a unscented transform x(k + 1) = f (x(k), u(k + 1))
%     kf_structure.X = f(kf_structure.X,...
%         [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz],...
%         G,T);
%     return;

x_previous = kf_structure.X;
P_previous = kf_structure.P;
u = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
Pu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
kf_structure.X = localization_filter_model_f('evaluate',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
df_dx = localization_filter_model_f('df_dx',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
df_du = localization_filter_model_f('df_du',x_previous,u,G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
kf_structure.P = df_dx*P_previous*df_dx' + df_du*Pu*df_du' + kf_structure.Q_ekf_decoupled;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [kf_structure] = localization_filter_prediction_ukf2(kf_structure, imumeasure, G, T)
%Etapa de predi��o usando a unscented transform x(k + 1) = f (x(k), u(k + 1))
%     kf_structure.X = f(kf_structure.X,...
%         [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz],...
%         G,T);
%     return;

Xm = [kf_structure.X; imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
Pu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
Px = [kf_structure.P zeros(kf_structure.Nstates,6); zeros(6,kf_structure.Nstates) Pu];
DimX = kf_structure.Nstates + 6;
DimY = kf_structure.Nstates;
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
   YSamples(:,i)      = localization_filter_model_f('evaluate',XSamples(1:kf_structure.Nstates,i),XSamples((kf_structure.Nstates+1):(kf_structure.Nstates+6),i),G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
   YSamples(:,i+DimX) = localization_filter_model_f('evaluate',XSamples(1:kf_structure.Nstates,i+DimX),XSamples((kf_structure.Nstates+1):(kf_structure.Nstates+6),i+DimX),G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias);
end,
WSamples(1,UtNSamples) = kappa/(DimX+kappa);
XSamples(:,UtNSamples) = Xm;
YSamples(:,UtNSamples) = localization_filter_model_f('evaluate',XSamples(1:kf_structure.Nstates,UtNSamples),XSamples((kf_structure.Nstates+1):(kf_structure.Nstates+6),UtNSamples),G,T,kf_structure.flagestimateaccelerometerbias, kf_structure.flagestimategyrometerbias); 

Ym 	= zeros(DimY,1);
Py  = zeros(DimY,DimY);
for ns=1:UtNSamples
   Ym = Ym + WSamples(1,ns)*YSamples(:,ns);
end,
for ns=1:UtNSamples
   Py = Py + WSamples(1,ns)*(YSamples(:,ns)-Ym)*(YSamples(:,ns)-Ym)';
end,

kf_structure.X(1:4) = quaternions_correctsign(Ym(1:4), kf_structure.X(1:4));
kf_structure.Xsigma = YSamples; % guarda os sigma pontos do estado para serem aproveitados na corre�ao.
kf_structure.P = Py + kf_structure.Q_ukf;


