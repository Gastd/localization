%Função que simula a aquisição de dados por uma IMU
function [imumeasure] = imu(vehiclestate,G,flagnoise)

[imumeasure] = imumeasure_init(flagnoise);

if isstruct(vehiclestate)==0
    return;
end

%Dado os angulos pela trajetória, calcular a velocidade angular
imumeasure.wx = vehiclestate.droll_dt - vehiclestate.dyaw_dt*sin(vehiclestate.pitch);
imumeasure.wy = vehiclestate.dpitch_dt*cos(vehiclestate.roll) + vehiclestate.dyaw_dt*cos(vehiclestate.pitch)*sin(vehiclestate.roll);
imumeasure.wz = - vehiclestate.dpitch_dt*sin(vehiclestate.roll) + vehiclestate.dyaw_dt*cos(vehiclestate.pitch)*cos(vehiclestate.roll);

%Desvio Padrão do gyros de 0.0085º/s em x, 0.0075º/s em y,
%0.0213º/s em z
%adicionando os erros de medição
imumeasure.wx = imumeasure.wx + flagnoise*sqrt(imumeasure.wxvariance)*randn(1);
imumeasure.wy = imumeasure.wy + flagnoise*sqrt(imumeasure.wyvariance)*randn(1);
imumeasure.wz = imumeasure.wz + flagnoise*sqrt(imumeasure.wzvariance)*randn(1);

%Dado a matriz de rotação e possível calcular no sistema do corpo qual a
%aceleração incluindo a gravidade
Reb = euler2dcm([vehiclestate.roll,vehiclestate.pitch,vehiclestate.yaw]');
a = Reb'*([vehiclestate.d2x_dt2; vehiclestate.d2y_dt2; vehiclestate.d2z_dt2] - G); % gravidade no eixo z

imumeasure.ax = a(1);
imumeasure.ay = a(2);
imumeasure.az = a(3);

%desvio Padrão nos acelerômetros de 0.0659 m/s^2 em x, 0.0639 m/s^2 em y,
%0.0213 m/s^2 em z
%adicionando os erros de medição
imumeasure.ax = imumeasure.ax + flagnoise*sqrt(imumeasure.axvariance)*randn(1);
imumeasure.ay = imumeasure.ay + flagnoise*sqrt(imumeasure.ayvariance)*randn(1);
imumeasure.az = imumeasure.az + flagnoise*sqrt(imumeasure.azvariance)*randn(1);
