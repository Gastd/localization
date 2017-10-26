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
	LocalizationFilter filter_struct;
	//First step: initialize the filter structure
	localization_init(LOCALIZATION_ALGORITHMCODE_EKF_DECOUPLED, 1, &filter_struct);

	//Second step: initialize the covariance matrices and initial state estimate
	PGMATRIX_DATA(filter_struct.pX, X_q0, 1) = 1.0;
	PGMATRIX_DATA(filter_struct.pX, X_q1, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_q2, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_q3, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_px, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_py, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_pz, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_vx, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_vy, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_vz, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_bax, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_bay, 1) = 0.0;
	PGMATRIX_DATA(filter_struct.pX, X_baz, 1) = 0.0;

	PGMATRIX_IDENTITY(filter_struct.pP);
	PGMATRIX_MULTIPLY_CONST(filter_struct.pP, 1000.0);

	PGMATRIX_IDENTITY(filter_struct.pPreset);
	PGMATRIX_MULTIPLY_CONST(filter_struct.pPreset, 1000.0);

	//Now, create the necessary variables
	GMATRIX_DECLARE(G, 3, 1); //Gravity vector
	GMATRIX_DATA(G, 1, 1) = 0.0;
	GMATRIX_DATA(G, 2, 1) = 0.0;
	GMATRIX_DATA(G, 3, 1) = 9.8182;

	GMATRIX_DECLARE(M, 3, 1); //Magnetic field vector
	GMATRIX_DATA(M, 1, 1) = 20.80224;
	GMATRIX_DATA(M, 2, 1) = -7.8082;
	GMATRIX_DATA(M, 3, 1) = -8.63213;

	SonarMeasure sonar_measure;
	sonar_init(&sonar_measure);

	ImuMeasure imu_measure;
	imu_init(&imu_measure);

	GpsMeasure gps_measure;
	gps_init(&gps_measure);

	MagnetometerMeasure magnetometer_measure;
	magnetometer_init(&magnetometer_measure);

	double T=0.1;

	//And put some data into them
	//We don't use the sonar for now
	//Fake measurements for IMU for prediction
	imu_measure.ax = 1.0;
	imu_measure.ay = 0.0;
	imu_measure.az = -10.0;
	imu_measure.wx = 0.1;
	imu_measure.wy = -0.05;
	imu_measure.wz = 0.02;

	localization_filter_prediction(&filter_struct, &imu_measure, &G, T);

	//Fake GPS and Mag data
	gps_measure.FlagValidPositionMeasure = 1;
	gps_measure.FlagValidVelocityMeasure = 1;
	PGMATRIX_DATA(gps_measure.pPosition, 1, 1) = 1.0;
	PGMATRIX_DATA(gps_measure.pPosition, 2, 1) = 1.0;
	PGMATRIX_DATA(gps_measure.pPosition, 3, 1) = 1.0;
	PGMATRIX_DATA(gps_measure.pVelocity, 1, 1) = 0.4;
	PGMATRIX_DATA(gps_measure.pVelocity, 2, 1) = 0.3;
	PGMATRIX_DATA(gps_measure.pVelocity, 3, 1) = 0.7;

	magnetometer_measure.FlagValidMeasure = 1;
	magnetometer_measure.mx = 26.22;
	magnetometer_measure.my = 3.081;
	magnetometer_measure.mz = -18.73;
	localization_filter_correction(&filter_struct, &gps_measure, &imu_measure, &magnetometer_measure, &sonar_measure, &M, &G, T);

	PGMATRIX_PRINT_MATLABFORM(filter_struct.pX);

	magnetometer_close(&magnetometer_measure);
	gps_close(&gps_measure);
	imu_close(&imu_measure);
	sonar_close(&sonar_measure);
	localization_close(&filter_struct);

	return 1;
}
