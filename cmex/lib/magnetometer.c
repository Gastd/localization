#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "magnetometer.h"

int magnetometer_init(PMAGNETOMETERMEASURE pMagnetometerMeasure)
{
	pMagnetometerMeasure->FlagValidMeasure = 0;

	pMagnetometerMeasure->mxvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_MAGNETOMETER_X_STANDART_DEVIATION);
	pMagnetometerMeasure->myvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_MAGNETOMETER_Y_STANDART_DEVIATION);
	pMagnetometerMeasure->mzvariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_MAGNETOMETER_Z_STANDART_DEVIATION);

    return 1;
}                      

int magnetometer_close(PMAGNETOMETERMEASURE pMagnetometerMeasure)
{
    return 1;
}
