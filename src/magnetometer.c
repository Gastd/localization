#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "magnetometer.h"

int magnetometer_init(MagnetometerMeasure *magnetometer_measure_ptr)
{
	magnetometer_measure_ptr->FlagValidMeasure = 0;

	magnetometer_measure_ptr->mxvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_MAGNETOMETER_X_STANDART_DEVIATION);
	magnetometer_measure_ptr->myvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_MAGNETOMETER_Y_STANDART_DEVIATION);
	magnetometer_measure_ptr->mzvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_MAGNETOMETER_Z_STANDART_DEVIATION);

    return 1;
}                      

int magnetometer_close(MagnetometerMeasure *magnetometer_measure_ptr)
{
    return 1;
}
