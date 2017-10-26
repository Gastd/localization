#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "gps.h"

int gps_init(GpsMeasure *gps_measure_ptr)
{
	gps_measure_ptr->FlagValidPositionMeasure = 0;
	gps_measure_ptr->FlagValidVelocityMeasure = 0;
	gps_measure_ptr->pPosition = PGMATRIX_ALLOC(3,1);
	gps_measure_ptr->pVelocity = PGMATRIX_ALLOC(3,1);
	gps_measure_ptr->pPPosition = PGMATRIX_ALLOC(3,3);
	gps_measure_ptr->pPVelocity = PGMATRIX_ALLOC(3,3);

	PGMATRIX_ZEROES(gps_measure_ptr->pPPosition);
	PGMATRIX_DATA(gps_measure_ptr->pPPosition,1,1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_X_STANDART_DEVIATION);
	PGMATRIX_DATA(gps_measure_ptr->pPPosition,2,2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_Y_STANDART_DEVIATION);
	PGMATRIX_DATA(gps_measure_ptr->pPPosition,3,3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_Z_STANDART_DEVIATION);

	PGMATRIX_ZEROES(gps_measure_ptr->pPVelocity);
	PGMATRIX_DATA(gps_measure_ptr->pPVelocity,1,1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_VX_STANDART_DEVIATION);
	PGMATRIX_DATA(gps_measure_ptr->pPVelocity,2,2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_VY_STANDART_DEVIATION);
	PGMATRIX_DATA(gps_measure_ptr->pPVelocity,3,3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_VZ_STANDART_DEVIATION);

    return 1;
}                      

int gps_close(GpsMeasure *gps_measure_ptr)
{
	PGMATRIX_FREE(gps_measure_ptr->pPosition);
	PGMATRIX_FREE(gps_measure_ptr->pVelocity);
	PGMATRIX_FREE(gps_measure_ptr->pPPosition);
	PGMATRIX_FREE(gps_measure_ptr->pPVelocity);

    return 1;
}                      
