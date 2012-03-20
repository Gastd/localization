/**************************************************************
 Programa: mexquaternions2euler.c
 Objetivo: Versao Cmex de quaternions2euler
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
	PQUATERNIONS pq;

	/* Verificar as entradas */
	if(nrhs!=1){
		mexErrMsgTxt("One input required.");
		}
	if(nlhs!=1){
		mexErrMsgTxt("One output required.");
		}

	/* Variaveis de entrada */
    pq = PGMATRIX_ALLOC_FROM_MXARRAY(prhs[0]);

    /* Função principal */
    rotation_quaternions2euler(pq, &Roll, &Pitch, &Yaw);
    
    /* Variaveis de saida */
	plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);	pData = mxGetPr(plhs[0]);
    pData[0] = Roll;
    pData[1] = Pitch;
    pData[2] = Yaw;
}
	
