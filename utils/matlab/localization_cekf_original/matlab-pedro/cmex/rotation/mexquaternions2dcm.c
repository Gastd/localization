/**************************************************************
 Programa: mexquaternions2dcm.c
 Objetivo: Versao Cmex de quaternions2dcm
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
	DCM_DECLARE(R);
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
    rotation_quaternions2dcm(pq, &R);
    
    /* Variaveis de saida */
    plhs[0] = PGMATRIX_COPY_TO_MXARRAY(&R);
}
	
