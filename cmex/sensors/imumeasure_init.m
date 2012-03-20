%Fun��o que simula a aquisi��o de dados por uma IMU
function [imumeasure] = imumeasure_init(flagnoise)

imumeasure.ax = 0;
imumeasure.ay = 0;
imumeasure.az = 0;
imumeasure.wx = 0;
imumeasure.wy = 0;
imumeasure.wz = 0;

variance_wx = 100*((0.0085)*pi/180)^2;
variance_wy = 100*((0.0075)*pi/180)^2;
variance_wz = 100*((0.0213)*pi/180)^2;
variance_ax = 100*((0.0659))^2;
variance_ay = 100*((0.0639))^2;
variance_az = 100*((0.0624))^2;

imumeasure.axvariance = flagnoise*variance_ax;
imumeasure.ayvariance = flagnoise*variance_ay;
imumeasure.azvariance = flagnoise*variance_az;
imumeasure.wxvariance = flagnoise*variance_wx;
imumeasure.wyvariance = flagnoise*variance_wy;
imumeasure.wzvariance = flagnoise*variance_wz;

