/**********************************************************************************************
***********************************************************************************************
***  File:			GMatrix_matlab.c
***	 Author:		Geovany Araujo Borges
***	 Contents:		Implementation of Matrix type and related functions for use in C and C++.
***					Matlab interfacing functions
***********************************************************************************************
**********************************************************************************************/
//#include "stdafx.h" // MSVC6.0 may require this.
#include "mex.h"

#include "gmatrix.h"
#include "gmatrix_matlab.h"

/**********************************************************************************************
***** GMatrix: Includes.
**********************************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/**********************************************************************************************
***** GMatrix: CMEX
**********************************************************************************************/

PGMATRIX PGMATRIX_ALLOC_FROM_MXARRAY(const mxArray *pmxarray)
{
	PGMATRIX mat;
	double *pData;
	int i,j;

	pData = mxGetPr(pmxarray);
	mat = PGMATRIX_ALLOC(mxGetM(pmxarray),mxGetN(pmxarray));
	for(i=1;i<=mat->Nr;++i){
		for(j=1;j<=mat->Nc;++j){
			PGMATRIX_DATA(mat,i,j) = pData[(i-1)+(j-1)*mat->Nr];
		}
	}

	return(mat);
}

mxArray *PGMATRIX_COPY_TO_MXARRAY(PGMATRIX pMat)
{
	mxArray *pmxarray;
	double *pData;
	int i,j;

	pmxarray = mxCreateDoubleMatrix(pMat->Nr,pMat->Nc,mxREAL);

	pData = mxGetPr(pmxarray);
	for(i=1;i<=pMat->Nr;++i){
		for(j=1;j<=pMat->Nc;++j){
			pData[(i-1)+(j-1)*pMat->Nr] = PGMATRIX_DATA(pMat,i,j);
		}
	}

	return(pmxarray);
}

