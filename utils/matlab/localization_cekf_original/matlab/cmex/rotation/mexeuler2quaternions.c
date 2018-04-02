/**************************************************************
 Programa: mexeuler2quaternions.c
 Objetivo: Versao Cmex de euler2quaternions
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
	QUATERNIONS_DECLARE(q);

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
    rotation_euler2quaternions(Roll, Pitch, Yaw, &q);
    
    /* Variaveis de saida */
	plhs[0] = mxCreateDoubleMatrix(4,1,mxREAL);	pData = mxGetPr(plhs[0]);
    pData[0] = QUATERNIONS_Q0(q);
    pData[1] = QUATERNIONS_Q1(q);
    pData[2] = QUATERNIONS_Q2(q);
    pData[3] = QUATERNIONS_Q3(q);
}
	
