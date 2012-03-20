/**************************************************************
 Programa: mexdcm2euler.c
 Objetivo: Versao Cmex de dcm2euler
**************************************************************/
/* A ne pas oublier: x[i][j] (Format C) <=> x[i+j*M] (Format CMEX) */

#include "mex.h"
#include <math.h>

#include "gmatrix.h"
#include "gmatrix_matlab.h"
#include "gmatrix_linalg.h"
#include "rotation.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int		i,j;
	double  Roll, Pitch, Yaw, *pData;
	PDCM pR;

	/* Verificar as entradas */
	if(nrhs!=1){
		mexErrMsgTxt("One input required.");
		}
	if(nlhs!=1){
		mexErrMsgTxt("One output required.");
		}

	if ( (mxGetM(prhs[0])!=mxGetN(prhs[0])) ){
		mexErrMsgTxt("R must be square.");
	}
	if ( (mxGetM(prhs[0])!= 3) ){
		mexErrMsgTxt("R must be 3 x 3.");
	}
	/* Variaveis de entrada */
    pR = PGMATRIX_ALLOC_FROM_MXARRAY(prhs[0]);

    /* Função principal */
//    *pw = CI_Optimization(gamma, &Pa, &Pb, &DummyMatrices);
    rotation_dcm2euler(pR, &Roll, &Pitch, &Yaw);
    
    /* Variaveis de saida */
	plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);	pData = mxGetPr(plhs[0]);
    pData[0] = Roll;
    pData[1] = Pitch;
    pData[2] = Yaw;

    PDCM_FREE(pR);
}
	
