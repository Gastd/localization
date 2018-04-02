/**************************************************************
 Programa: mexquaternions_correctsign.c
 Objetivo: Versao Cmex de quaternions_correctsign
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
	PQUATERNIONS pqcurrent;
	PQUATERNIONS pqprevious;

	/* Verificar as entradas */
	if(nrhs!=2){
		mexErrMsgTxt("Two inputs required.");
		}
	if(nlhs!=1){
		mexErrMsgTxt("One output required.");
		}

	/* Variaveis de entrada */
    pqcurrent  = PGMATRIX_ALLOC_FROM_MXARRAY(prhs[0]);
    pqprevious = PGMATRIX_ALLOC_FROM_MXARRAY(prhs[1]);

    /* Função principal */
    rotation_quaternionscorrectsign(pqcurrent, pqprevious);
    
    /* Variaveis de saida */
	plhs[0] = PGMATRIX_COPY_TO_MXARRAY(pqcurrent);
}
	
