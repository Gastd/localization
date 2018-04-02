%fun��o que calcula a atitude apartir dos dados dos magnet�metros e do
%aceler�metro
function q = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted)

FlagConsiderPrediction = 0;
if nargin > 4
    FlagConsiderPrediction = 1;
end

mag = [magnetometermeasure.mx; magnetometermeasure.my; magnetometermeasure.mz];

%armazena a acelera��o (for�a espec�fica). A acerela��o � negada, pois leva-se em conta que os
%acelerometros medem a for�a espec�fica sob sua massa e n�o a acelera��o
a = zeros(3,1);
a(1) = -imumeasure.ax;
a(2) = -imumeasure.ay;
a(3) = -imumeasure.az;

%normaliza��o dos vetores de acelera��o e de campo magn�tico
mag = (mag/sqrt(mag(1)^2 + mag(2)^2 + mag(3)^2));
a = a/sqrt(a(1)^2 + a(2)^2 + a(3)^2);
G = G/sqrt(G(1)^2 + G(2)^2 + G(3)^2);
M = M/sqrt(M(1)^2 + M(2)^2 + M(3)^2);


% Inicializa as variaveis do corpo. sendo "a" as medidas do acelerometro e
% "mag" as medidas do magnetometro. Lembrando que a fun��o cross(,)
% representa produto vetorial entre dois vetores.
aux = a + mag;
I_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

aux = cross(I_b,(a - mag));
J_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

K_b = cross(I_b,J_b);

%Inicializa as variaveis do sistema inercial
aux = G + M;
I_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

aux = cross(I_i,(G - M));
J_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);

K_i = cross(I_i,J_i);

% Calculo da matriz de rota��o pelo m�todo TRIAD melhorado
R_i_b = ([I_b, J_b, K_b]*([I_i, J_i, K_i]'));

%Cria��o do quat�rnio de rota��o. (Deve se observar que a matriz de
%rota��o � do sistema de refer�ncia para o sistema do corpo a contr�rio do
%que o Padilha escreveu. Dados para a confer�ncia est�o em
%http://www.uel.br/proppg/semina/pdf/semina_28_1_22_19.pdf).
R_b_i = R_i_b';
q = dcm2quaternions(R_b_i);

if FlagConsiderPrediction,
    q = quaternions_correctsign(q, q_predicted);
end