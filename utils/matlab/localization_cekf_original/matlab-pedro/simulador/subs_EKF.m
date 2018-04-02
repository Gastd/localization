%%% Funcao auxiliar que apenas atribui valores numericos aos paramentros do
%%% EKF.


function [H,h_est] = subs_EKF(kf_structure,M)

%Atribuicao de valor as variaveis simbolicas
mx  = M(1,1);
my  = M(2,1);
mz  = M(3,1);

q0  =  kf_structure.X(1,1);
q1  =  kf_structure.X(2,1);
q2  =  kf_structure.X(3,1);
q3  =  kf_structure.X(4,1);

h1 = (q0^2+q1^2-q2^2-q3^2)*mx + 2*(q1*q2+q0*q3)*my + 2*(q1*q3-q0*q2)*mz;
h2 = 2*(q1*q2-q0*q3)*mx + (q0^2-q1^2+q2^2-q3^2)*my + 2*(q2*q3+q0*q1)*mz;
h3 = 2*(q1*q3+q0*q2)*mx + 2*(q2*q3-q0*q1)*my + (q0^2-q1^2-q2^2+q3^2)*mz;

%Estimativa de leitura baseada na predicao do estado e na funcao
%nao-linear h_mag
h_est = [h1 h2 h3]';

% Jacobiana da funcao de medicao avaliada na ultima estimativa
H = zeros(3,13);

H(1,1) = 2*( q0*mx + q3*my - q2*mz );
H(1,2) = 2*( q1*mx + q2*my + q3*mz );
H(1,3) = 2*( -q2*mx + q1*my - q0*mz );
H(1,4) = 2*( -q3*mx + q0*my + q1*mz );

H(2,1) = 2*( -q3*mx + q0*my + q1*mz );
H(2,2) = 2*( q2*mx - q1*my + q0*mz );
H(2,3) = 2*( q1*mx + q2*my + q3*mz );
H(2,4) = 2*( -q0*mx - q3*my + q2*mz );

H(3,1) = 2*( q2*mx - q1*my + q0*mz );
H(3,2) = 2*( q3*mx - q0*my - q1*mz );
H(3,3) = 2*( q0*mx + q3*my - q2*mz );
H(3,4) = 2*( q1*mx + q2*my + q3*mz );

     

