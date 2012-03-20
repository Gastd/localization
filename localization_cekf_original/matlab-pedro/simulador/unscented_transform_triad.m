%Transformada Unscented

%   Parametros de entrada
%       Xm -> vetor em torno do qual são calculados os pontos sigma
%       (medidas do acelerometro e do magnetometro)
%       Px -> matriz de covariancia associada a Xm (variancias das medicoes
%       experimentais).
%
%   Parametros de saida
%       Ym -> media transformada
%       Py -> matriz de covariancia transformada

function [Ym,Py] = unscented_transform_triad(imumeasure,magnetometermeasure,M,G,q_predicted)

%%% Propagação de incertezas por Unscented Transform.

%Vetor de medicoes
Xm = [imumeasure.ax imumeasure.ay imumeasure.az magnetometermeasure.mx magnetometermeasure.my magnetometermeasure.mz]';

%Acelerometro
Paccel = [imumeasure.axvariance 0 0; 0 imumeasure.ayvariance 0 ; 0 0 imumeasure.azvariance];
%Magnetometro
Pmag = [magnetometermeasure.mxvariance 0 0; 0 magnetometermeasure.myvariance 0 ; 0 0 magnetometermeasure.mzvariance];
%Matriz de covariancia das medidas
Px = [Paccel zeros(3); zeros(3) Pmag];

%Dimensao do vetor de entradas (3 medidas do acelerometro e 3 do
%magnetometro)
DimX = 6;
% Dimensao do quaternio de saida
DimY = 4;

%Parametro misterioso do Prof. Geovany
kappa = 1;

%Numero de amostras para a transformada unscented
UtNSamples = 2*DimX+1;

%Estrutura para armazenar os pontos sigma
XSamples = zeros(DimX,UtNSamples);
%Estrutura para armazenar os p ontos sigma transformados
YSamples = zeros(DimY,UtNSamples);
% Estrutura para armazenar os pesos dos pontos sigma
WSamples = zeros(1,UtNSamples);

%Vetor para geracao dos pontos sigma
Q = sqrtm((DimX+kappa)*Px)';

for i=1:DimX
    
    %Peso dos pontos sigma i e i+DimX
    WSamples(1,i)   = 1/(2*(DimX+kappa));
    WSamples(1,i+DimX) = 1/(2*(DimX+kappa));
    
    %Ponto sigma i
    XSamples(:,i)   = Xm + Q(:,i);
    
    %Ponto sigma i+DimX
    XSamples(:,i+DimX) = Xm - Q(:,i);
        
    %Armazena os pontos sigma para serem transformados
    imumeasure.ax           = XSamples(1,i);
    imumeasure.ay           = XSamples(2,i);
    imumeasure.az           = XSamples(3,i);
    magnetometermeasure.mx  = XSamples(4,i);
    magnetometermeasure.my  = XSamples(5,i);
    magnetometermeasure.mz  = XSamples(6,i);   
    YSamples(:,i) = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    
    imumeasure.ax           = XSamples(1,i+DimX);
    imumeasure.ay           = XSamples(2,i+DimX);
    imumeasure.az           = XSamples(3,i+DimX);
    magnetometermeasure.mx  = XSamples(4,i+DimX);
    magnetometermeasure.my  = XSamples(5,i+DimX);
    magnetometermeasure.mz  = XSamples(6,i+DimX);   
    YSamples(:,i+DimX) = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
end

%Transformacao do ultimo ponto, Xm
WSamples(1,UtNSamples) = kappa/(DimX+kappa);
XSamples(:,UtNSamples) = Xm;

imumeasure.ax           = Xm(1,1);
imumeasure.ay           = Xm(2,1);
imumeasure.az           = Xm(3,1);
magnetometermeasure.mx  = Xm(4,1);
magnetometermeasure.my  = Xm(5,1);
magnetometermeasure.mz  = Xm(6,1);
YSamples(:,UtNSamples) = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);

%Calculo do quaternio transformado e de sua matriz de covariancias
Ym 	= zeros(DimY,1);
Py  = zeros(DimY,DimY);

%Quaternio transformado, cujo calculo eh feito por meio da ponderacao das
%transformacoes dos pontos sigma.
for ns=1:UtNSamples
   Ym = Ym + WSamples(1,ns)*YSamples(:,ns);
end

%Matriz de covariancias da transformacao
for ns=1:UtNSamples
    dy = YSamples(:,ns)-Ym;    
    Py = Py + WSamples(1,ns)*(dy)*(dy)';
end
