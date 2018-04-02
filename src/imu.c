#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "imu.h"

int imu_init(ImuMeasure *imu_measure_ptr)
{
	imu_measure_ptr->FlagValidAccerometerMeasure = 0;
	imu_measure_ptr->FlagValidGyrometerMeasure = 0;

	imu_measure_ptr->axvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_ACCELEROMETER_X_STANDART_DEVIATION);
	imu_measure_ptr->ayvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_ACCELEROMETER_Y_STANDART_DEVIATION);
	imu_measure_ptr->azvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_ACCELEROMETER_Z_STANDART_DEVIATION);

	imu_measure_ptr->wxvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GYROMETER_X_STANDART_DEVIATION);
	imu_measure_ptr->wyvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GYROMETER_Y_STANDART_DEVIATION);
	imu_measure_ptr->wzvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GYROMETER_Z_STANDART_DEVIATION);

	return 1;
}

int imu_close(ImuMeasure *imu_measure_ptr)
{
	return 1;
}
