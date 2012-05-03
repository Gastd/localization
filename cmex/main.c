#include "lib/gmatrix.h"
#include "lib/rotation.h"
#include "lib/gps.h"
#include "lib/imu.h"
#include "lib/kalman.h"
#include "lib/magnetometer.h"
#include "lib/sonar.h"
#include "lib/localization.h"

int main(void)
{
	LOCALIZATIONFILTERSTRUCT FilterStruct;
	//First step: initialize the filter structure
	localization_init(LOCALIZATION_ALGORITHMCODE_EKF_DECOUPLED, 1, &FilterStruct);

	//Second step: initialize the covariance matrices and initial state estimate
	PGMATRIX_DATA(FilterStruct.pX, X_q0, 1) = 1.0;
	PGMATRIX_DATA(FilterStruct.pX, X_q1, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_q2, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_q3, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_px, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_py, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_pz, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_vx, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_vy, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_vz, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_bax, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_bay, 1) = 0.0;
	PGMATRIX_DATA(FilterStruct.pX, X_baz, 1) = 0.0;

	PGMATRIX_IDENTITY(FilterStruct.pP);
	PGMATRIX_MULTIPLY_CONST(FilterStruct.pP, 1000.0);

	PGMATRIX_IDENTITY(FilterStruct.pPreset);
	PGMATRIX_MULTIPLY_CONST(FilterStruct.pPreset, 1000.0);

	//Now, create the necessary variables
	GMATRIX_DECLARE(G, 3, 1); //Gravity vector
	GMATRIX_DATA(G, 1, 1) = 0.0;
	GMATRIX_DATA(G, 2, 1) = 0.0;
	GMATRIX_DATA(G, 3, 1) = 9.8182;

	GMATRIX_DECLARE(M, 3, 1); //Magnetic field vector
	GMATRIX_DATA(M, 1, 1) = 20.80224;
	GMATRIX_DATA(M, 2, 1) = -7.8082;
	GMATRIX_DATA(M, 3, 1) = -8.63213;

	SONARMEASURE SonarMeasure;
	sonar_init(&SonarMeasure);

	IMUMEASURE IMUMeasure;
	imu_init(&IMUMeasure);

	GPSMEASURE GPSMeasure;
	gps_init(&GPSMeasure);

	MAGNETOMETERMEASURE MagnetometerMeasure;
	magnetometer_init(&MagnetometerMeasure);

	double T=0.1;

	//And put some data into them
	//We don't use the sonar for now
	//Fake measurements for IMU for prediction
	IMUMeasure.ax = 1.0;
	IMUMeasure.ay = 0.0;
	IMUMeasure.az = -10.0;
	IMUMeasure.wx = 0.1;
	IMUMeasure.wy = -0.05;
	IMUMeasure.wz = 0.02;

	localization_filter_prediction(&FilterStruct, &IMUMeasure, &G, T);

	//Fake GPS and Mag data
	GPSMeasure.FlagValidPositionMeasure = 1;
	GPSMeasure.FlagValidVelocityMeasure = 1;
	PGMATRIX_DATA(GPSMeasure.pPosition, 1, 1) = 1.0;
	PGMATRIX_DATA(GPSMeasure.pPosition, 2, 1) = 1.0;
	PGMATRIX_DATA(GPSMeasure.pPosition, 3, 1) = 1.0;
	PGMATRIX_DATA(GPSMeasure.pVelocity, 1, 1) = 0.4;
	PGMATRIX_DATA(GPSMeasure.pVelocity, 2, 1) = 0.3;
	PGMATRIX_DATA(GPSMeasure.pVelocity, 3, 1) = 0.7;

	MagnetometerMeasure.FlagValidMeasure = 1;
	MagnetometerMeasure.mx = 26.22;
	MagnetometerMeasure.my = 3.081;
	MagnetometerMeasure.mz = -18.73;
	localization_filter_correction(&FilterStruct, &GPSMeasure, &IMUMeasure, &MagnetometerMeasure, &SonarMeasure, &M, &G, T);

	PGMATRIX_PRINT_MATLABFORM(FilterStruct.pX);

	magnetometer_close(&MagnetometerMeasure);
	gps_close(&GPSMeasure);
	imu_close(&IMUMeasure);
	sonar_close(&SonarMeasure);
	localization_close(&FilterStruct);

	return 1;
}
