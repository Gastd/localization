#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "imu.h"

int imu_init(PIMUMEASURE pIMUMeasure)
{
	pIMUMeasure->FlagValidAccerometerMeasure = 0;
	pIMUMeasure->FlagValidGyrometerMeasure = 0;

	pIMUMeasure->axvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_ACCELEROMETER_X_STANDART_DEVIATION);
	pIMUMeasure->ayvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_ACCELEROMETER_Y_STANDART_DEVIATION);
	pIMUMeasure->azvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_ACCELEROMETER_Z_STANDART_DEVIATION);

	pIMUMeasure->wxvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GYROMETER_X_STANDART_DEVIATION);
	pIMUMeasure->wyvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GYROMETER_Y_STANDART_DEVIATION);
	pIMUMeasure->wzvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GYROMETER_Z_STANDART_DEVIATION);

	return 1;
}

int imu_close(PIMUMEASURE pIMUMeasure)
{
	return 1;
}
