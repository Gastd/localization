#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "gps.h"

int gps_init(PGPSMEASURE pGPSMeasure)
{
	pGPSMeasure->FlagValidPositionMeasure = 0;
	pGPSMeasure->FlagValidVelocityMeasure = 0;
	pGPSMeasure->pPosition = PGMATRIX_ALLOC(3,1);
	pGPSMeasure->pVelocity = PGMATRIX_ALLOC(3,1);
	pGPSMeasure->pPPosition = PGMATRIX_ALLOC(3,3);
	pGPSMeasure->pPVelocity = PGMATRIX_ALLOC(3,3);

	PGMATRIX_ZEROES(pGPSMeasure->pPPosition);
	PGMATRIX_DATA(pGPSMeasure->pPPosition,1,1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_X_STANDART_DEVIATION);
	PGMATRIX_DATA(pGPSMeasure->pPPosition,2,2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_Y_STANDART_DEVIATION);
	PGMATRIX_DATA(pGPSMeasure->pPPosition,3,3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_Z_STANDART_DEVIATION);

	PGMATRIX_ZEROES(pGPSMeasure->pPVelocity);
	PGMATRIX_DATA(pGPSMeasure->pPVelocity,1,1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_VX_STANDART_DEVIATION);
	PGMATRIX_DATA(pGPSMeasure->pPVelocity,2,2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_VY_STANDART_DEVIATION);
	PGMATRIX_DATA(pGPSMeasure->pPVelocity,3,3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_GPS_VZ_STANDART_DEVIATION);

    return 1;
}                      

int gps_close(PGPSMEASURE pGPSMeasure)
{
	PGMATRIX_FREE(pGPSMeasure->pPosition);
	PGMATRIX_FREE(pGPSMeasure->pVelocity);
	PGMATRIX_FREE(pGPSMeasure->pPPosition);
	PGMATRIX_FREE(pGPSMeasure->pPVelocity);

    return 1;
}                      
