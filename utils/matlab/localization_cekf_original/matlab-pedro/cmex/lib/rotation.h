/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: rotation.h 
Conteudo: Cabe�alho de fun��es em c�digo C de manipula��o de matrizes de cosseno diretores (DCM).
Autor: G. A. Borges.
Atualiza��es: 
	- 28/07/2008: cria��o do exemplo, por Geovany A. Borges
*****************************************************************************/

#ifndef ROTATION_H
#define ROTATION_H

// Definicoes de uso externo
#define PQUATERNIONS				PGMATRIX

#define QUATERNIONS_DECLARE(q)		GMATRIX_DECLARE(q,4,1)

#define QUATERNIONS_Q0(q)			PQUATERNIONS_Q0((&q))
#define QUATERNIONS_Q1(q)			PQUATERNIONS_Q1((&q))
#define QUATERNIONS_Q2(q)			PQUATERNIONS_Q2((&q))
#define QUATERNIONS_Q3(q)			PQUATERNIONS_Q3((&q))

#define PQUATERNIONS_Q0(pq)			PGMATRIX_DATA(pq,1,1)
#define PQUATERNIONS_Q1(pq)			PGMATRIX_DATA(pq,2,1)
#define PQUATERNIONS_Q2(pq)			PGMATRIX_DATA(pq,3,1)
#define PQUATERNIONS_Q3(pq)			PGMATRIX_DATA(pq,4,1)

#define PDCM			PGMATRIX
#define PDCM_ALLOC(pR)	pR = PGMATRIX_ALLOC(3,3)
#define PDCM_FREE(pR)	PGMATRIX_FREE(pR)
#define DCM_DECLARE(R)	GMATRIX_DECLARE(R,3,3)

#define DCM_R11(R)			GMATRIX_DATA(R,1,1)
#define DCM_R12(R)			GMATRIX_DATA(R,1,2)
#define DCM_R13(R)			GMATRIX_DATA(R,1,3)
#define DCM_R21(R)			GMATRIX_DATA(R,2,1)
#define DCM_R22(R)			GMATRIX_DATA(R,2,2)
#define DCM_R23(R)			GMATRIX_DATA(R,2,3)
#define DCM_R31(R)			GMATRIX_DATA(R,3,1)
#define DCM_R32(R)			GMATRIX_DATA(R,3,2)
#define DCM_R33(R)			GMATRIX_DATA(R,3,3)

#define PDCM_R11(pR)		PGMATRIX_DATA(pR,1,1)
#define PDCM_R12(pR)		PGMATRIX_DATA(pR,1,2)
#define PDCM_R13(pR)		PGMATRIX_DATA(pR,1,3)
#define PDCM_R21(pR)		PGMATRIX_DATA(pR,2,1)
#define PDCM_R22(pR)		PGMATRIX_DATA(pR,2,2)
#define PDCM_R23(pR)		PGMATRIX_DATA(pR,2,3)
#define PDCM_R31(pR)		PGMATRIX_DATA(pR,3,1)
#define PDCM_R32(pR)		PGMATRIX_DATA(pR,3,2)
#define PDCM_R33(pR)		PGMATRIX_DATA(pR,3,3)


// Prototipos externos:
int rotation_dcm2euler(PDCM pR, double *proll, double *ppitch, double *pyaw);
int rotation_euler2dcm(double roll, double pitch, double yaw, PDCM pR);
int rotation_quaternionsnormalize(PQUATERNIONS pq);
int rotation_quaternions2dcm(PQUATERNIONS pq, PDCM pR);
int rotation_quaternions2euler(PQUATERNIONS pq, double *proll, double *ppitch, double *pyaw);
int rotation_dcm2quaternions(PDCM pR, PQUATERNIONS pq);
int rotation_quaternionscorrectsign(PQUATERNIONS pqcurrent, PQUATERNIONS pqprevious);

#endif
