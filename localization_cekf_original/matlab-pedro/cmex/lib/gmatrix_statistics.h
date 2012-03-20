/**********************************************************************************************
***********************************************************************************************
***  File:			GMatrix_statistics.h
***	 Author:		Geovany Araujo Borges
***	 Contents:		GMatrix_statistics header file.
***********************************************************************************************
**********************************************************************************************/


/**********************************************************************************************
***** GMatrix: Function prototypes
**********************************************************************************************/

#define GMATRIX_MEAN(Mat) PGMATRIX_MEAN(&Mat)
double PGMATRIX_MEAN(PGMATRIX pMat);

#define GMATRIX_VARIANCE(Mat) PGMATRIX_VARIANCE(&Mat)
double PGMATRIX_VARIANCE(PGMATRIX pMat);

#define GMATRIX_COVARIANCE_PROPAGATION_COPY(MatCovResult, MatGain, MatCov, MatDummy) PGMATRIX_COVARIANCE_PROPAGATION_COPY(&MatCovResult, &MatGain, &MatCov, &MatDummy)
void PGMATRIX_COVARIANCE_PROPAGATION_COPY(PGMATRIX pMatCovResult, PGMATRIX pMatGain, PGMATRIX pMatCov, PGMATRIX pMatDummy);

#define GMATRIX_COVARIANCE_PROPAGATION_ADD(MatCovResult, MatGain, MatCov, MatDummy) PGMATRIX_COVARIANCE_PROPAGATION_ADD(&MatCovResult, &MatGain, &MatCov, &MatDummy)
void PGMATRIX_COVARIANCE_PROPAGATION_ADD(PGMATRIX pMatCovResult, PGMATRIX pMatGain, PGMATRIX pMatCov, PGMATRIX pMatDummy);

#define GMATRIX_MAHALANOBIS_DISTANCE(MatResidual, MatCovResidual, MatDummy1, MatDummy2) PGMATRIX_MAHALANOBIS_DISTANCE(&MatResidual, &MatCovResidual, &MatDummy1, &MatDummy2)
double PGMATRIX_MAHALANOBIS_DISTANCE(PGMATRIX pMatResidual, PGMATRIX pMatCovResidual, PGMATRIX pMatDummy1, PGMATRIX pMatDummy2);

