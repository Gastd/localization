#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "filter_parameters.h"
#include "sonar.h"

int sonar_init(SonarMeasure *sonar_measure_ptr)
{
	sonar_measure_ptr->FlagValidMeasure = 0;
	sonar_measure_ptr->range = 0;
	sonar_measure_ptr->rangevariance = GMATRIXMACRO_SQR(FILTER_PARAMETERS_SONAR_STANDART_DEVIATION);

	sonar_measure_ptr->pR_s2b = PGMATRIX_ALLOC(3,3);
	sonar_measure_ptr->pt_s2b = PGMATRIX_ALLOC(3,1);

	PGMATRIX_IDENTITY(sonar_measure_ptr->pR_s2b);
	PGMATRIX_ZEROES(sonar_measure_ptr->pt_s2b);

    return 1; 
}                      

int sonar_close(SonarMeasure *sonar_measure_ptr)
{
	PGMATRIX_FREE(sonar_measure_ptr->pR_s2b);
	PGMATRIX_FREE(sonar_measure_ptr->pt_s2b);

    return 1; 
}                      

