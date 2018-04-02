#ifndef UNSCENTEDTRANSFORM_H
#define UNSCENTEDTRANSFORM_H

//#include "mex.h"
// Requires GMatrix.
#include <math.h>
#include <stdlib.h>

/**********************************************************************************************
***** Uncertainty propagation: USCENTED TRANSFORMATION
**********************************************************************************************/
#define UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM_INIT(DimY,DimX)\
	GMATRIX_DECLARE(XMean		,(DimX)	,1);\
	GMATRIX_DECLARE(YMean		,(DimY)	,1);\
	GMATRIX_DECLARE(XCov			,(DimX)	,(DimX));\
	GMATRIX_DECLARE(YCov			,(DimY)	,(DimY));\
	GMATRIX_DECLARE(_utX			,(DimX)	,1);\
	GMATRIX_DECLARE(_utY			,(DimY)	,1);\
	GMATRIX_DECLARE(_utXSamples	,(DimX)	,(2*(DimX)+1));\
	GMATRIX_DECLARE(_utYSamples	,(DimY)	,(2*(DimX)+1));\
	GMATRIX_DECLARE(_utWSamples	,(1)	,(2*(DimX)+1));\
	GMATRIX_DECLARE(_utESamples	,(DimY)	,1);\
	GMATRIX_DECLARE(_utQSamples	,(DimX)	,(DimX))

#define UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM(FunctionCall) _UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM(&XMean,&YMean,&XCov,&YCov,&_utX,&_utY,&_utXSamples,&_utYSamples,&_utWSamples,&_utESamples,&_utQSamples,&Parameters,FunctionCall)
	
__inline void _UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM(
	PGMATRIX pXMean, 
	PGMATRIX pYMean, 
	PGMATRIX pXCov, 
	PGMATRIX pYCov, 
	PGMATRIX pX, 
	PGMATRIX pY, 
	PGMATRIX pXSamples, 
	PGMATRIX pYSamples, 
	PGMATRIX pWSamples, 
	PGMATRIX pESamples, 
	PGMATRIX pQSamples,
	PGMATRIX pParameters,
	void pFunctionCall(PGMATRIX pY, PGMATRIX pX, PGMATRIX pParameters) )
{
	int	ns,i,j,DimX,DimY,NSamples;
	double kappa;

	// Init parameters:
	DimX = pXMean->Nr;
	DimY = pYMean->Nr;
	NSamples = 2*DimX+1;
	kappa = 1.0;

	// Init QSamples = cholesky((DimX+kappa)*XCov)
	PGMATRIX_MULTIPLY_CONST(pXCov,(DimX+kappa)); // XCov is no more valid!
	PGMATRIX_CHOLESKY(pQSamples, pXCov);
	// Generate samples:
	for (i=1;i<=DimX;++i){
		PGMATRIX_DATA(pWSamples,1,i)		= 1/(2.0*(DimX+kappa));
		PGMATRIX_DATA(pWSamples,1,i+DimX)	= 1/(2.0*(DimX+kappa));
		for (j=1;j<=DimX;++j){
			PGMATRIX_DATA(pXSamples,j,i)		= PGMATRIX_DATA(pXMean,j,1) + PGMATRIX_DATA(pQSamples,j,i);
			PGMATRIX_DATA(pXSamples,j,i+DimX)	= PGMATRIX_DATA(pXMean,j,1) - PGMATRIX_DATA(pQSamples,j,i);
		}
	}
	PGMATRIX_DATA(pWSamples,1,NSamples) = kappa/(double)(DimX+kappa);
	PGMATRIX_SUBMATRIX_COPY(pXSamples,1,NSamples,pXMean);
	for (i=1;i<=NSamples;++i){
		PGMATRIX_COPY_COLUMN(pX, 1, pXSamples, i);
		pFunctionCall(pY,pX,pParameters);
		PGMATRIX_COPY_COLUMN(pYSamples, i, pY, 1);
	}
	
	// Compute estimates:
	PGMATRIX_ZEROES(pYMean);
	PGMATRIX_ZEROES(pYCov);
	for (ns=1;ns<=NSamples;++ns){
		for (i=1;i<=DimY;++i){
			PGMATRIX_DATA(pYMean,i,1) += PGMATRIX_DATA(pWSamples,1,ns)*PGMATRIX_DATA(pYSamples,i,ns);
		}
	}
	for (ns=1;ns<=NSamples;++ns){
		for (i=1;i<=DimY;++i){
			PGMATRIX_DATA(pESamples,i,1) = PGMATRIX_DATA(pYSamples,i,ns) - PGMATRIX_DATA(pYMean,i,1);
		}
		for (i=1;i<=DimY;++i){
			for (j=1;j<=DimY;++j){
				PGMATRIX_DATA(pYCov,i,j) += PGMATRIX_DATA(pWSamples,1,ns)*PGMATRIX_DATA(pESamples,i,1)*PGMATRIX_DATA(pESamples,j,1);
			}
		}
	}
}

#endif //UNSCENTEDTRANSFORM_H
