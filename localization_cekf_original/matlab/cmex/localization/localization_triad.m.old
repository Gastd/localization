%função que calcula a atitude apartir dos dados dos magnetômetros e do
%acelerômetro
function q = estimator_triad(imumeasure, magnetometermeasure, M, G, q_predicted)

FlagConsiderPrediction = 0;
if nargin > 4
    FlagConsiderPrediction = 1;
end

mag = [magnetometermeasure.mx; magnetometermeasure.my; magnetometermeasure.mz];

%armazena a aceleração (força específica). A acerelação é negada, pois leva-se em conta que os
%acelerometros medem a força específica sob sua massa e não a aceleração
a = zeros(3,1);
a(1) = -imumeasure.ax;
a(2) = -imumeasure.ay;
a(3) = -imumeasure.az;

%normalização dos vetores de aceleração e de campo magnético
mag = (mag/sqrt(mag(1)^2 + mag(2)^2 + mag(3)^2));
a = a/sqrt(a(1)^2 + a(2)^2 + a(3)^2);
G = G/sqrt(G(1)^2 + G(2)^2 + G(3)^2);
M = M/sqrt(M(1)^2 + M(2)^2 + M(3)^2);


% Inicializa as variaveis do corpo. sendo "a" as medidas do acelerometro e
% "mag" as medidas do magnetometro. Lembrando que a função cross(,)
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

% Calculo da matriz de rotação pelo método TRIAD melhorado
R_i_b = ([I_b, J_b, K_b]*([I_i, J_i, K_i]'));

%Criação do quatérnio de rotação. (Deve se observar que a matriz de
%rotação é do sistema de referência para o sistema do corpo a contrário do
%que o Padilha escreveu. Dados para a conferência estão em
%http://www.uel.br/proppg/semina/pdf/semina_28_1_22_19.pdf).
R_b_i = R_i_b';
q = dcm2quaternions(R_b_i);

if FlagConsiderPrediction,
    q = quaternions_correctsign(q, q_predicted);
end