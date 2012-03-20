#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "sonar.h"

int sonar_init(PSONARMEASURE pSonarMeasure)
{
	pSonarMeasure->FlagValidMeasure = 0;
	pSonarMeasure->range = 0;
	pSonarMeasure->rangevariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_SONAR_STANDART_DEVIATION);

	pSonarMeasure->pR_s2b = PGMATRIX_ALLOC(3,3);
	pSonarMeasure->pt_s2b = PGMATRIX_ALLOC(3,1);

	PGMATRIX_IDENTITY(pSonarMeasure->pR_s2b);
	PGMATRIX_ZEROES(pSonarMeasure->pt_s2b);

    return 1; 
}                      

int sonar_close(PSONARMEASURE pSonarMeasure)
{
	PGMATRIX_FREE(pSonarMeasure->pR_s2b);
	PGMATRIX_FREE(pSonarMeasure->pt_s2b);

    return 1; 
}                      

