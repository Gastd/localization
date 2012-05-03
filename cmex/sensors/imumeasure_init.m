%Fun��o que simula a aquisi��o de dados por uma IMU
function [imumeasure] = imumeasure_init(flagnoise)

imumeasure.ax = 0;
imumeasure.ay = 0;
imumeasure.az = 0;
imumeasure.wx = 0;
imumeasure.wy = 0;
imumeasure.wz = 0;

variance_wx = 0.0001;
variance_wy = 0.0001;
variance_wz = 0.0001;
variance_ax = 0.01;
variance_ay = 0.01;
variance_az = 0.01;

imumeasure.axvariance = flagnoise*variance_ax;
imumeasure.ayvariance = flagnoise*variance_ay;
imumeasure.azvariance = flagnoise*variance_az;
imumeasure.wxvariance = flagnoise*variance_wx;
imumeasure.wyvariance = flagnoise*variance_wy;
imumeasure.wzvariance = flagnoise*variance_wz;

