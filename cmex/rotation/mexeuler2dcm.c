/**************************************************************
 Programa: mexeuler2dcm.c
 Objetivo: Versao Cmex de euler2dcm
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
	DCM_DECLARE(R);

	/* Verificar as entradas */
	if(nrhs!=1){
		mexErrMsgTxt("One input required.");
		}
	if(nlhs!=1){
		mexErrMsgTxt("One output required.");
		}

	/* Variaveis de entrada */
    pData  = mxGetPr(prhs[0]);
    Roll = pData[0];
    Pitch = pData[1];
    Yaw = pData[2];

    /* Função principal */
    rotation_euler2dcm(Roll, Pitch, Yaw, &R);
    
    /* Variaveis de saida */
    plhs[0] = PGMATRIX_COPY_TO_MXARRAY(&R);
}
	
