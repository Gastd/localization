#ifndef FILTER_PARAMETERS_H
#define FILTER_PARAMETERS_H

//GPS parameters (meters)
//#define FILTER_PARAMETERS_GPS_X_STANDART_DEVIATION						2.0
//#define FILTER_PARAMETERS_GPS_Y_STANDART_DEVIATION						2.0
//#define FILTER_PARAMETERS_GPS_Z_STANDART_DEVIATION						3.0
#define FILTER_PARAMETERS_GPS_X_STANDART_DEVIATION                      10.0
#define FILTER_PARAMETERS_GPS_Y_STANDART_DEVIATION                      10.0
#define FILTER_PARAMETERS_GPS_Z_STANDART_DEVIATION                      10.0
#define FILTER_PARAMETERS_GPS_VX_STANDART_DEVIATION						0.1
#define FILTER_PARAMETERS_GPS_VY_STANDART_DEVIATION						0.1
#define FILTER_PARAMETERS_GPS_VZ_STANDART_DEVIATION						0.2

//IMU parameters (meters, radians)
#define FILTER_PARAMETERS_ACCELEROMETER_X_STANDART_DEVIATION			0.1
#define FILTER_PARAMETERS_ACCELEROMETER_Y_STANDART_DEVIATION			0.1
#define FILTER_PARAMETERS_ACCELEROMETER_Z_STANDART_DEVIATION			0.1
#define FILTER_PARAMETERS_GYROMETER_X_STANDART_DEVIATION				10.0*(0.0085)*M_PI/180.0
#define FILTER_PARAMETERS_GYROMETER_Y_STANDART_DEVIATION				10.0*(0.0075)*M_PI/180.0
#define FILTER_PARAMETERS_GYROMETER_Z_STANDART_DEVIATION				10.0*(0.0213)*M_PI/180.0

//Magnetometer parameters (uTesla)
#define FILTER_PARAMETERS_MAGNETOMETER_X_STANDART_DEVIATION				0.5
#define FILTER_PARAMETERS_MAGNETOMETER_Y_STANDART_DEVIATION				0.5
#define FILTER_PARAMETERS_MAGNETOMETER_Z_STANDART_DEVIATION				0.5

//Sonar parameters
#define FILTER_PARAMETERS_SONAR_STANDART_DEVIATION						0.01

//Q Matrix
//Note: this WILL BE SQUARED! That is, this is the standard deviation, not the variance
#define FILTER_PARAMETERS_Q_MATRIX_Q0									0.01
#define FILTER_PARAMETERS_Q_MATRIX_Q1									0.01
#define FILTER_PARAMETERS_Q_MATRIX_Q2									0.01
#define FILTER_PARAMETERS_Q_MATRIX_Q3									0.01
#define FILTER_PARAMETERS_Q_MATRIX_POSITION_X							0.01
#define FILTER_PARAMETERS_Q_MATRIX_POSITION_Y							0.01
#define FILTER_PARAMETERS_Q_MATRIX_POSITION_Z							0.01
#define FILTER_PARAMETERS_Q_MATRIX_VELOCITY_X							0.0001
#define FILTER_PARAMETERS_Q_MATRIX_VELOCITY_Y							0.0001
#define FILTER_PARAMETERS_Q_MATRIX_VELOCITY_Z							0.0001
#define FILTER_PARAMETERS_Q_MATRIX_ACCELEROMETER_BIAS_X					0.001
#define FILTER_PARAMETERS_Q_MATRIX_ACCELEROMETER_BIAS_Y					0.001
#define FILTER_PARAMETERS_Q_MATRIX_ACCELEROMETER_BIAS_Z					0.001

//R - pseudomeasurementnorm
#define FILTER_PARAMETERS_R_PSEUDOMEASUREMENTNORM						1e-5

//R - convertedmeasurementtriad
#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_1_1			0.001
#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_2_2			0.001
#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_3_3			0.001
#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_4_4			0.001

#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_1_1				0.01
#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_2_2				0.01
#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_3_3				0.01
#define FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_4_4				0.01

#define FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_1_1             2.0
#define FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_2_2             2.0
#define FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_3_3             2.0

#define FILTER_PARAMETERS_R_MAGNETOMETER_EKF_DECOUPLED_1_1              0.1
#define FILTER_PARAMETERS_R_MAGNETOMETER_EKF_DECOUPLED_2_2              0.1
#define FILTER_PARAMETERS_R_MAGNETOMETER_EKF_DECOUPLED_3_3              0.1

#endif //FILTER_PARAMETERS_H
