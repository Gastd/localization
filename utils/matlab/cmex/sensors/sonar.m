% Fun��o que simula o sonar
function sonarmeasure = sonar(vehiclestate, flagnoise)

sonarmeasure = sonarmeasure_init(flagnoise);

if isstruct(vehiclestate)==0
    return;
end

%Cria a matriz de rota��o dado o quat�rnio
R = euler2dcm([vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw]');

%Armazena o desvio padr�o de 0,01 m para o ru�do

%calcula os valores das supostas medidas do altimetro
sonarmeasure.range = vehiclestate.z / R(3,3);     %a dist�ncia medida pelo altimetro deve ser igual a altura do corpo(z no sistema de referencia inercial)
                                               %divido pelo elemento a33 da
                                               %matriz de rotacao
%adiciona o ruido caso habilitado
sonarmeasure.range = sonarmeasure.range + flagnoise*sqrt(sonarmeasure.rangevariance)*randn(1);