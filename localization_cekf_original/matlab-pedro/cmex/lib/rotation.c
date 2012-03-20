/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: rotation.c 
Conteudo: Funções em código C de manipulação de matrizes de cosseno diretores (DCM).
Autor: G. A. Borges.
Atualizações: 
	- 19/09/2008: criação por Geovany A. Borges, modulo sem init/close
*****************************************************************************/
// Cabecalhos de bibliotecas C run-time
#include <math.h>
#include <stdio.h>

// Cabecalhos de modulos do projeto
#include "gmatrix.h"
#include "gmatrix_linalg.h"
#include "rotation.h"

// Definicoes de uso interno

// Prototipos de funcoes internas ao modulo

// Variaveis globais do módulo

/*****************************************************************************
******************************************************************************
** FUNCOES COM CHAMADA EXTERNA
******************************************************************************
*****************************************************************************/

/*****************************************************************************
*** int rotation_quaternionscorrectsign(PQUATERNIONS pqcurrent, PQUATERNIONS pqprevious)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_quaternionscorrectsign(PQUATERNIONS pqcurrent, PQUATERNIONS pqprevious)
{
	QUATERNIONS_DECLARE(plus);
	QUATERNIONS_DECLARE(minus);
	
	// now we check if we want q or -q
	// plus  = q_predicted - q;
	PGMATRIX_SUBTRACT_COPY(&plus, pqprevious, pqcurrent);
	// minus = q_predicted + q;
	PGMATRIX_ADD_COPY(&minus, pqprevious, pqcurrent);
	// if (norm(plus) > norm(minus))
	//		q = - q;
	// end
	if( GMATRIX_NORM(plus) > GMATRIX_NORM(minus) ){
		PGMATRIX_MULTIPLY_CONST(pqcurrent, -1.0);
	}
	
	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int rotation_quaternions2euler(PQUATERNIONS pq, double *proll, double *ppitch, double *pyaw)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_quaternions2euler(PQUATERNIONS pq, double *proll, double *ppitch, double *pyaw)
{
	double norm;
	DCM_DECLARE(R);

	// make quaternion unit length: done in quaternions2dcm
	//	R = quaternions2dcm(q);
	if(!rotation_quaternions2dcm(pq, &R)) return 0;

	//	[roll,pitch,yaw] = dcm2euler(R);
	if(!rotation_dcm2euler(&R, proll, ppitch, pyaw)) return 0;

	
	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int rotation_euler2quaternions(PDCM pR, double proll, double ppitch, double pyaw)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_euler2quaternions(double roll, double pitch, double yaw, PQUATERNIONS pq)
{
	DCM_DECLARE(R);

	// R = euler2dcm(roll,pitch,yaw);
	if(!rotation_euler2dcm(roll, pitch, yaw, &R)) return 0;

	// q = dcm2quaternions(R);
	if(!rotation_dcm2quaternions(&R, pq)) return 0;
	
	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int rotation_quaternions2dcm(PQUATERNIONS pq, PDCM pR)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_quaternions2dcm(PQUATERNIONS pq, PDCM pR)
{
	double a;

	// make quaternion unit length:
	rotation_quaternionsnormalize(pq);

	// in acordance with Titterton & Weston, pg. 48, and Coutsias & Romero, last pg
	// R = 2*[(q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2)/2    (q(2)*q(3)-q(1)*q(4))                   (q(2)*q(4)+q(1)*q(3));
	PDCM_R11(pR) = 1.0*(PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq) + PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq) - PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq) - PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R12(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q2(pq) - PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R13(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q3(pq) + PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q2(pq));

	//	   (q(2)*q(3)+q(1)*q(4))                    (q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2)/2   (q(3)*q(4)-q(1)*q(2));
	PDCM_R21(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q2(pq) + PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R22(pR) = 1.0*(PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq) - PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq) + PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq) - PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R23(pR) = 2.0*(PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q3(pq) - PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q1(pq));

	//	   (q(2)*q(4)-q(1)*q(3))                    (q(3)*q(4)+q(1)*q(2))                   (q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2)/2];
	PDCM_R31(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q3(pq) - PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q2(pq));
	PDCM_R32(pR) = 2.0*(PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q3(pq) + PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q1(pq));
	PDCM_R33(pR) = 1.0*(PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq) - PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq) - PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq) + PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq));

	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int rotation_dcm2quaternions(PDCM pR, PQUATERNIONS pq)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_dcm2quaternions(PDCM pR, PQUATERNIONS pq)
{
	double a;
//	% improved code: Antonio Padilha
//	q0 = .5*sqrt(max(0,1+R(1,1)+R(2,2)+R(3,3)));
	a = 1.0 + PGMATRIX_DATA(pR,1,1) + PGMATRIX_DATA(pR,2,2) + PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q0(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q0(pq) = 0.0;

//	q1 = .5*sqrt(max(0,1+R(1,1)-R(2,2)-R(3,3)));
	a = 1.0 + PGMATRIX_DATA(pR,1,1) - PGMATRIX_DATA(pR,2,2) - PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q1(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q1(pq) = 0.0;

//	q2 = .5*sqrt(max(0,1-R(1,1)+R(2,2)-R(3,3)));
	a = 1.0 - PGMATRIX_DATA(pR,1,1) + PGMATRIX_DATA(pR,2,2) - PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q2(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q2(pq) = 0.0;

//	q3 = .5*sqrt(max(0,1-R(1,1)-R(2,2)+R(3,3)));
	a = 1.0 - PGMATRIX_DATA(pR,1,1) - PGMATRIX_DATA(pR,2,2) + PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q3(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q3(pq) = 0.0;

//	q1 = abs(q1)*sign(R(3,2)-R(2,3));
	if ((PGMATRIX_DATA(pR,3,2)-PGMATRIX_DATA(pR,2,3)) >= 0.0)
		PQUATERNIONS_Q1(pq) =  fabs(PQUATERNIONS_Q1(pq));
	else
		PQUATERNIONS_Q1(pq) = -fabs(PQUATERNIONS_Q1(pq));

//	q2 = abs(q2)*sign(R(1,3)-R(3,1));
	if ((PGMATRIX_DATA(pR,1,3)-PGMATRIX_DATA(pR,3,1)) >= 0.0)
		PQUATERNIONS_Q2(pq) =  fabs(PQUATERNIONS_Q2(pq));
	else
		PQUATERNIONS_Q2(pq) = -fabs(PQUATERNIONS_Q2(pq));

//	q3 = abs(q3)*sign(R(2,1)-R(1,2));
	if ((PGMATRIX_DATA(pR,2,1)-PGMATRIX_DATA(pR,1,2)) >= 0.0)
		PQUATERNIONS_Q3(pq) =  fabs(PQUATERNIONS_Q3(pq));
	else
		PQUATERNIONS_Q3(pq) = -fabs(PQUATERNIONS_Q3(pq));

	// Retorna 
    return 1; 
}             

/*****************************************************************************
*** int rotation_quaternions2euler(PQUATERNIONS pq, double *proll, double *ppitch, double *pyaw)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_quaternionsnormalize(PQUATERNIONS pq)
{
	double norm;
	DCM_DECLARE(R);

	// % make quaternion unit length:
	// q = q / norm(q);
	norm  = PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq);
	norm += PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq);
	norm += PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq);
	norm += PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq);
	norm =  sqrt(norm);

	if(norm == 0.0) return 0;

	PQUATERNIONS_Q0(pq) = PQUATERNIONS_Q0(pq) / norm;
	PQUATERNIONS_Q1(pq) = PQUATERNIONS_Q1(pq) / norm;
	PQUATERNIONS_Q2(pq) = PQUATERNIONS_Q2(pq) / norm;
	PQUATERNIONS_Q3(pq) = PQUATERNIONS_Q3(pq) / norm;
	
	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int rotation_dcm2euler(PDCM pR, double proll, double ppitch, double pyaw)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_dcm2euler(PDCM pR, double *proll, double *ppitch, double *pyaw)
{
	double cosPitch;

	// pitch = asin(-R(3,1));
	*ppitch = asin(-PGMATRIX_DATA(pR,3,1));
	
	cosPitch = cos(*ppitch);
	if (cosPitch == 0)
		return 0;
	
	//roll = atan2(R(3,2)/cos(pitch),R(3,3)/cos(pitch));
	*proll  = atan2(PGMATRIX_DATA(pR,3,2)/cosPitch,PGMATRIX_DATA(pR,3,3)/cosPitch);
	
	// yaw = atan2(R(2,1)/cos(pitch),R(1,1)/cos(pitch));
	*pyaw  = atan2(PGMATRIX_DATA(pR,2,1)/cosPitch,PGMATRIX_DATA(pR,1,1)/cosPitch);
	
	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int rotation_euler2dcm(double roll, double pitch, double yaw, PDCM pR)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int rotation_euler2dcm(double roll, double pitch, double yaw, PDCM pR)
{

	if (pR->Nr !=3){
		return 0;
	}
	if (pR->Nc !=3){
		return 0;
	}

	// R(1,:) = [cos(pitch)*cos(yaw) -cos(roll)*sin(yaw)+sin(roll)*sin(pitch)*cos(yaw)  sin(roll)*sin(yaw)+cos(roll)*sin(pitch)*cos(yaw)];
	PGMATRIX_DATA(pR,1,1) = cos(pitch)*cos(yaw);
	PGMATRIX_DATA(pR,1,2) = -cos(roll)*sin(yaw)+sin(roll)*sin(pitch)*cos(yaw);
	PGMATRIX_DATA(pR,1,3) = sin(roll)*sin(yaw)+cos(roll)*sin(pitch)*cos(yaw);

	// R(2,:) = [cos(pitch)*sin(yaw)  cos(roll)*cos(yaw)+sin(roll)*sin(pitch)*sin(yaw) -sin(roll)*cos(yaw)+cos(roll)*sin(pitch)*sin(yaw)];
	PGMATRIX_DATA(pR,2,1) = cos(pitch)*sin(yaw);
	PGMATRIX_DATA(pR,2,2) = cos(roll)*cos(yaw)+sin(roll)*sin(pitch)*sin(yaw); 
	PGMATRIX_DATA(pR,2,3) = -sin(roll)*cos(yaw)+cos(roll)*sin(pitch)*sin(yaw);
	
	// R(3,:) = [-sin(pitch)          sin(roll)*cos(pitch)                              cos(roll)*cos(pitch)];
	PGMATRIX_DATA(pR,3,1) = -sin(pitch);
	PGMATRIX_DATA(pR,3,2) = sin(roll)*cos(pitch);
	PGMATRIX_DATA(pR,3,3) = cos(roll)*cos(pitch);

	// Retorna 
    return 1; 
}                      

/*
void vecCross (PGMATRIX pResult, PGMATRIX pP1, PGMATRIX pP2)
{
	vecSkew(pga33_1_1,pP1);
	PGMATRIX_MULTIPLY_COPY (pResult,pga33_1_1,pP2);
}

void quatSet (PGMATRIX pQ,double q0,double q1,double q2,double q3)
{
	PGMATRIX_DATA(pQ,1,1) = q0;
	PGMATRIX_DATA(pQ,2,1) = q1;
	PGMATRIX_DATA(pQ,3,1) = q2;
	PGMATRIX_DATA(pQ,4,1) = q3;
}

float quatNorm(PGMATRIX pQ)
{
	return ( sqrt(PGMATRIX_DATA(pQ,1,1)*PGMATRIX_DATA(pQ,1,1) + PGMATRIX_DATA(pQ,2,1)*PGMATRIX_DATA(pQ,2,1) + PGMATRIX_DATA(pQ,3,1)*PGMATRIX_DATA(pQ,3,1) + PGMATRIX_DATA(pQ,4,1)*PGMATRIX_DATA(pQ,4,1)) );
}

void quatNormalize(PGMATRIX pQ)
{
	float norm = quatNorm(pQ);

	PGMATRIX_DATA(pQ,1,1) = PGMATRIX_DATA(pQ,1,1)/norm;
	PGMATRIX_DATA(pQ,2,1) = PGMATRIX_DATA(pQ,2,1)/norm;
	PGMATRIX_DATA(pQ,3,1) = PGMATRIX_DATA(pQ,3,1)/norm;
	PGMATRIX_DATA(pQ,4,1) = PGMATRIX_DATA(pQ,4,1)/norm;
}

void quatSkew (PGMATRIX pResult, float wx, float wy, float wz)
{
	// merwe's version
	PGMATRIX_DATA(pResult,1,1) = 0;
	PGMATRIX_DATA(pResult,1,2) = wx;
	PGMATRIX_DATA(pResult,1,3) = wy;
	PGMATRIX_DATA(pResult,1,4) = wz;
	PGMATRIX_DATA(pResult,2,1) = -wx;
	PGMATRIX_DATA(pResult,2,2) = 0;
	PGMATRIX_DATA(pResult,2,3) = -wz;
	PGMATRIX_DATA(pResult,2,4) = wy;
	PGMATRIX_DATA(pResult,3,1) = -wy;
	PGMATRIX_DATA(pResult,3,2) = wz;
	PGMATRIX_DATA(pResult,3,3) = 0;
	PGMATRIX_DATA(pResult,3,4) = -wx;
	PGMATRIX_DATA(pResult,4,1) = -wz;
	PGMATRIX_DATA(pResult,4,2) = -wy;
	PGMATRIX_DATA(pResult,4,3) = wx;
	PGMATRIX_DATA(pResult,4,4) = 0;
}

void quatMultiplication (PGMATRIX pR, PGMATRIX pP, PGMATRIX pQ)
{
	PGMATRIX_DATA(pR,1,1) = PGMATRIX_DATA(pP,1,1)*PGMATRIX_DATA(pQ,1,1)-PGMATRIX_DATA(pP,2,1)*PGMATRIX_DATA(pQ,2,1) -PGMATRIX_DATA(pP,3,1)*PGMATRIX_DATA(pQ,3,1)-PGMATRIX_DATA(pP,4,1)*PGMATRIX_DATA(pQ,4,1);
	PGMATRIX_DATA(pR,2,1) = PGMATRIX_DATA(pP,2,1)*PGMATRIX_DATA(pQ,1,1)+PGMATRIX_DATA(pP,1,1)*PGMATRIX_DATA(pQ,2,1) +PGMATRIX_DATA(pP,3,1)*PGMATRIX_DATA(pQ,4,1)-PGMATRIX_DATA(pP,4,1)*PGMATRIX_DATA(pQ,3,1);
	PGMATRIX_DATA(pR,3,1) = PGMATRIX_DATA(pP,1,1)*PGMATRIX_DATA(pQ,3,1)-PGMATRIX_DATA(pP,2,1)*PGMATRIX_DATA(pQ,4,1) +PGMATRIX_DATA(pP,3,1)*PGMATRIX_DATA(pQ,1,1)+PGMATRIX_DATA(pP,4,1)*PGMATRIX_DATA(pQ,2,1);
	PGMATRIX_DATA(pR,4,1) = PGMATRIX_DATA(pP,1,1)*PGMATRIX_DATA(pQ,4,1)+PGMATRIX_DATA(pP,2,1)*PGMATRIX_DATA(pQ,3,1) -PGMATRIX_DATA(pP,3,1)*PGMATRIX_DATA(pQ,2,1)+PGMATRIX_DATA(pP,4,1)*PGMATRIX_DATA(pQ,1,1);
}

void quatTransformation (PGMATRIX pPn, PGMATRIX pQbn, PGMATRIX pPb)
{
	// ver dissertacao, secao A.3.2
	
	PGMATRIX_DATA(pga33_1_1,1,1) = 2*PGMATRIX_DATA(pQbn,1,1)-1+pow(PGMATRIX_DATA(pQbn,2,1),2);
	PGMATRIX_DATA(pga33_1_1,1,2) = 2*PGMATRIX_DATA(pQbn,2,1)*PGMATRIX_DATA(pQbn,3,1) + 2*PGMATRIX_DATA(pQbn,1,1)*PGMATRIX_DATA(pQbn,4,1);
	PGMATRIX_DATA(pga33_1_1,1,3) = 2*PGMATRIX_DATA(pQbn,2,1)*PGMATRIX_DATA(pQbn,4,1) - 2*PGMATRIX_DATA(pQbn,1,1)*PGMATRIX_DATA(pQbn,3,1);
	PGMATRIX_DATA(pga33_1_1,2,1) = 2*PGMATRIX_DATA(pQbn,2,1)*PGMATRIX_DATA(pQbn,3,1) - 2*PGMATRIX_DATA(pQbn,1,1)*PGMATRIX_DATA(pQbn,4,1);
	PGMATRIX_DATA(pga33_1_1,2,2) = 2*PGMATRIX_DATA(pQbn,1,1)-1+pow(PGMATRIX_DATA(pQbn,3,1),2);
	PGMATRIX_DATA(pga33_1_1,2,3) = 2*PGMATRIX_DATA(pQbn,3,1)*PGMATRIX_DATA(pQbn,4,1) + 2*PGMATRIX_DATA(pQbn,1,1)*PGMATRIX_DATA(pQbn,2,1);
	PGMATRIX_DATA(pga33_1_1,3,1) = 2*PGMATRIX_DATA(pQbn,2,1)*PGMATRIX_DATA(pQbn,4,1) + 2*PGMATRIX_DATA(pQbn,1,1)*PGMATRIX_DATA(pQbn,3,1);
	PGMATRIX_DATA(pga33_1_1,3,2) = 2*PGMATRIX_DATA(pQbn,3,1)*PGMATRIX_DATA(pQbn,4,1) - 2*PGMATRIX_DATA(pQbn,1,1)*PGMATRIX_DATA(pQbn,2,1);
	PGMATRIX_DATA(pga33_1_1,3,3) = 2*PGMATRIX_DATA(pQbn,1,1)-1+pow(PGMATRIX_DATA(pQbn,4,1),2);
	
	PGMATRIX_MULTIPLY_COPY (pPn,pga33_1_1,pPb);
}

void quatPropagation (PGMATRIX pQ, PGMATRIX pQant, float v, float sX, float sY, float sZ)
{
	// procedimento descrito em merwe phd thesis
// 	PGMATRIX pA, pB, pBaux1, pBaux2;
// 	pA = PGMATRIX_ALLOC(4,4);
// 	pB = PGMATRIX_ALLOC(4,4);
// 	pBaux1 = PGMATRIX_ALLOC(4,4);
// 	pBaux2 = PGMATRIX_ALLOC(4,4);

	// primeiramente, calcula-se a exponencial, dada por B
// 	quatSkew(pA,sX,sY,sZ);
// 	PGMATRIX_IDENTITY(pB);
// 	PGMATRIX_MULTIPLY_CONST_COPY(pBaux1,pB,cos(v/2));
// 	PGMATRIX_MULTIPLY_CONST_COPY(pBaux2,pA,-sin(v/2)/v);
// 	PGMATRIX_ADD_COPY(pB,pBaux1,pBaux2);

	// em seguida, multiplica-se por qk para encontrar qk+1
// 	PGMATRIX_MULTIPLY_COPY (pQ, pB, pQant);

	// utilizando a eq desenvolvida, nao precisa das matrizes
	float cv2 = cos(v/2);
	float sv2v = (sin(v/2)/v);
	float Qant1 = PGMATRIX_DATA(pQant,1,1);
	float Qant2 = PGMATRIX_DATA(pQant,2,1);
	float Qant3 = PGMATRIX_DATA(pQant,3,1);
	float Qant4 = PGMATRIX_DATA(pQant,4,1);
	PGMATRIX_DATA(pQ,1,1) = cv2*Qant1 - sv2v * (sX*Qant2 + sY*Qant3 + sZ*Qant4);
	PGMATRIX_DATA(pQ,2,1) = cv2*Qant2 - sv2v * (-sX*Qant1 - sZ*Qant3 + sY*Qant4);
	PGMATRIX_DATA(pQ,3,1) = cv2*Qant3 - sv2v * (-sY*Qant1 + sZ*Qant2 - sX*Qant4);
	PGMATRIX_DATA(pQ,4,1) = cv2*Qant4 - sv2v * (-sZ*Qant1 - sY*Qant2 + sX*Qant3);

// 	PGMATRIX_FREE(pBaux1);
// 	PGMATRIX_FREE(pBaux2);
// 	PGMATRIX_FREE(pB);
// 	PGMATRIX_FREE(pA);
}

void dcm2quat (PGMATRIX pQ, PGMATRIX pC)
{
	// 2 methods may be used:

	// as in books Titterton, Applied and Quaternions, with singularity
	# if 0
	PGMATRIX_DATA(pQ,1,1) = .5 * sqrt(1 + PGMATRIX_DATA(pC,1,1) + PGMATRIX_DATA(pC,2,2) + PGMATRIX_DATA(pC,3,3));
	PGMATRIX_DATA(pQ,2,1) = 1/(4 * PGMATRIX_DATA(pQ,1,1)) * (PGMATRIX_DATA(pC,3,2) - PGMATRIX_DATA(pC,2,3));
	PGMATRIX_DATA(pQ,3,1) = 1/(4 * PGMATRIX_DATA(pQ,1,1)) * (PGMATRIX_DATA(pC,1,3) - PGMATRIX_DATA(pC,3,1));
	PGMATRIX_DATA(pQ,4,1) = 1/(4 * PGMATRIX_DATA(pQ,1,1)) * (PGMATRIX_DATA(pC,2,1) - PGMATRIX_DATA(pC,1,2));
	#endif
	
	// as proposed by euclidean space website
	#if 1
	PGMATRIX_DATA(pQ,1,1) = .5 * sqrt(fmaxf(0,1 + PGMATRIX_DATA(pC,1,1) + PGMATRIX_DATA(pC,2,2) + PGMATRIX_DATA(pC,3,3)));
	PGMATRIX_DATA(pQ,2,1) = .5 * sqrt(fmaxf(0,1 + PGMATRIX_DATA(pC,1,1) - PGMATRIX_DATA(pC,2,2) - PGMATRIX_DATA(pC,3,3)));
	PGMATRIX_DATA(pQ,3,1) = .5 * sqrt(fmaxf(0,1 - PGMATRIX_DATA(pC,1,1) + PGMATRIX_DATA(pC,2,2) - PGMATRIX_DATA(pC,3,3)));
	PGMATRIX_DATA(pQ,4,1) = .5 * sqrt(fmaxf(0,1 - PGMATRIX_DATA(pC,1,1) - PGMATRIX_DATA(pC,2,2) + PGMATRIX_DATA(pC,3,3)));

	PGMATRIX_DATA(pQ,2,1) = copysignf(PGMATRIX_DATA(pQ,2,1), PGMATRIX_DATA(pC,3,2) - PGMATRIX_DATA(pC,2,3));
	PGMATRIX_DATA(pQ,3,1) = copysignf(PGMATRIX_DATA(pQ,3,1), PGMATRIX_DATA(pC,1,3) - PGMATRIX_DATA(pC,3,1));
	PGMATRIX_DATA(pQ,4,1) = copysignf(PGMATRIX_DATA(pQ,4,1), PGMATRIX_DATA(pC,2,1) - PGMATRIX_DATA(pC,1,2));
	#endif
}

void dcmTransformation (PGMATRIX pPn, PGMATRIX pCbn, PGMATRIX pPb)
{
	PGMATRIX_MULTIPLY_COPY (pPn,pCbn,pPb);
}

void dcmFromAcelMag (PGMATRIX pCbi, PGMATRIX pFi, PGMATRIX pMi, PGMATRIX pFb, PGMATRIX pMb)
{
	// 2 methods may be used:

	// the simple proposed by Applied.. p.280
	// NOT YET ADAPTED TO NEW AUXILIAR VARIABLES
	#if 0
	PGMATRIX pTmpI,pTmpB,pTmpCross,pTmpInverse;
	pTmpCross = PGMATRIX_ALLOC(3,1);
	pTmpInverse = PGMATRIX_ALLOC(3,3);
	pTmpI = PGMATRIX_ALLOC(3,3);
	pTmpB = PGMATRIX_ALLOC(3,3);

	// normaliza vetores atuais
	vecNormalize(pgFb);
	vecNormalize(pgMb);

	vecCross(pTmpCross,pFi,pMi);

	PGMATRIX_DATA(pTmpI,1,1) = PGMATRIX_DATA(pFi,1,1);
	PGMATRIX_DATA(pTmpI,1,2) = PGMATRIX_DATA(pMi,1,1);
	PGMATRIX_DATA(pTmpI,1,3) = PGMATRIX_DATA(pTmpCross,1,1);
	PGMATRIX_DATA(pTmpI,2,1) = PGMATRIX_DATA(pFi,2,1);
	PGMATRIX_DATA(pTmpI,2,2) = PGMATRIX_DATA(pMi,2,1);
	PGMATRIX_DATA(pTmpI,2,3) = PGMATRIX_DATA(pTmpCross,2,1);
	PGMATRIX_DATA(pTmpI,3,1) = PGMATRIX_DATA(pFi,3,1);
	PGMATRIX_DATA(pTmpI,3,2) = PGMATRIX_DATA(pMi,3,1);
	PGMATRIX_DATA(pTmpI,3,3) = PGMATRIX_DATA(pTmpCross,3,1);

	vecCross(pTmpCross,pFb,pMb);

	PGMATRIX_DATA(pTmpB,1,1) = PGMATRIX_DATA(pFb,1,1);
	PGMATRIX_DATA(pTmpB,1,2) = PGMATRIX_DATA(pMb,1,1);
	PGMATRIX_DATA(pTmpB,1,3) = PGMATRIX_DATA(pTmpCross,1,1);
	PGMATRIX_DATA(pTmpB,2,1) = PGMATRIX_DATA(pFb,2,1);
	PGMATRIX_DATA(pTmpB,2,2) = PGMATRIX_DATA(pMb,2,1);
	PGMATRIX_DATA(pTmpB,2,3) = PGMATRIX_DATA(pTmpCross,2,1);
	PGMATRIX_DATA(pTmpB,3,1) = PGMATRIX_DATA(pFb,3,1);
	PGMATRIX_DATA(pTmpB,3,2) = PGMATRIX_DATA(pMb,3,1);
	PGMATRIX_DATA(pTmpB,3,3) = PGMATRIX_DATA(pTmpCross,3,1);

	PGMATRIX_INVERSE3_COPY(pTmpInverse,pTmpB);
	PGMATRIX_MULTIPLY_COPY(pCbi,pTmpI,pTmpInverse);

	PGMATRIX_FREE(pTmpCross);
	PGMATRIX_FREE(pTmpInverse);
	PGMATRIX_FREE(pTmpB);
	PGMATRIX_FREE(pTmpI);
	#endif

	// the improved TRIAD algorithm proposed by Yong Li
	#if 1

	// normaliza vetores atuais
	vecNormalize(pFb);
	vecNormalize(pMb);

	// B = [i_B|j_B|k_B], i_B = f_B+m_B/|f_B+m_B|, j_B = i_Bx(f_B-m_B)/|i_Bx(f_B-m_B)|, k_B = i_B x j_B
	PGMATRIX_ADD_COPY(pga31_2_5, pFb, pMb);
	PGMATRIX_MULTIPLY_CONST_COPY(pga31_2_1, pga31_2_5, 1/vecNorm(pga31_2_5));
	PGMATRIX_SUBSTRACT_COPY(pga31_2_5, pFb, pMb);
	vecCross(pga31_2_4, pga31_2_1, pga31_2_5);
	PGMATRIX_MULTIPLY_CONST_COPY(pga31_2_2, pga31_2_4, 1/vecNorm(pga31_2_4));	
	vecCross(pga31_2_3, pga31_2_1, pga31_2_2);

	PGMATRIX_DATA(pga33_2_2,1,1) = PGMATRIX_DATA(pga31_2_1,1,1);
	PGMATRIX_DATA(pga33_2_2,1,2) = PGMATRIX_DATA(pga31_2_2,1,1);
	PGMATRIX_DATA(pga33_2_2,1,3) = PGMATRIX_DATA(pga31_2_3,1,1);
	PGMATRIX_DATA(pga33_2_2,2,1) = PGMATRIX_DATA(pga31_2_1,2,1);
	PGMATRIX_DATA(pga33_2_2,2,2) = PGMATRIX_DATA(pga31_2_2,2,1);
	PGMATRIX_DATA(pga33_2_2,2,3) = PGMATRIX_DATA(pga31_2_3,2,1);
	PGMATRIX_DATA(pga33_2_2,3,1) = PGMATRIX_DATA(pga31_2_1,3,1);
	PGMATRIX_DATA(pga33_2_2,3,2) = PGMATRIX_DATA(pga31_2_2,3,1);
	PGMATRIX_DATA(pga33_2_2,3,3) = PGMATRIX_DATA(pga31_2_3,3,1);

	// I = analogous to B
	PGMATRIX_ADD_COPY(pga31_2_5, pFi, pMi);
	PGMATRIX_MULTIPLY_CONST_COPY(pga31_2_1, pga31_2_5, 1/vecNorm(pga31_2_5));
	PGMATRIX_SUBSTRACT_COPY(pga31_2_5, pFi, pMi);
	vecCross(pga31_2_4, pga31_2_1, pga31_2_5);
	PGMATRIX_MULTIPLY_CONST_COPY(pga31_2_2, pga31_2_4, 1/vecNorm(pga31_2_4));	
	vecCross(pga31_2_3, pga31_2_1, pga31_2_2);

	PGMATRIX_DATA(pga33_2_1,1,1) = PGMATRIX_DATA(pga31_2_1,1,1);
	PGMATRIX_DATA(pga33_2_1,1,2) = PGMATRIX_DATA(pga31_2_2,1,1);
	PGMATRIX_DATA(pga33_2_1,1,3) = PGMATRIX_DATA(pga31_2_3,1,1);
	PGMATRIX_DATA(pga33_2_1,2,1) = PGMATRIX_DATA(pga31_2_1,2,1);
	PGMATRIX_DATA(pga33_2_1,2,2) = PGMATRIX_DATA(pga31_2_2,2,1);
	PGMATRIX_DATA(pga33_2_1,2,3) = PGMATRIX_DATA(pga31_2_3,2,1);
	PGMATRIX_DATA(pga33_2_1,3,1) = PGMATRIX_DATA(pga31_2_1,3,1);
	PGMATRIX_DATA(pga33_2_1,3,2) = PGMATRIX_DATA(pga31_2_2,3,1);
	PGMATRIX_DATA(pga33_2_1,3,3) = PGMATRIX_DATA(pga31_2_3,3,1);

	PGMATRIX_TRANSPOSE_COPY(pga33_2_3,pga33_2_2);
	PGMATRIX_MULTIPLY_COPY(pCbi,pga33_2_1,pga33_2_3);
	#endif
}

void attitudePrediction (PGMATRIX pQ, PGMATRIX pQant, PGMATRIX pP, PGMATRIX pPant, PGMATRIX pQa, float wx, float wy, float wz, float dt)
{
	float v, sX, sY, sZ;
	sX = wx*dt;
	sY = wy*dt;
	sZ = wz*dt;
	v = sqrt(sX*sX + sY*sY + sZ*sZ);

	quatPropagation (pQ, pQant, v, sX, sY, sZ);

	// pga44_2_1 eh a jacobiana
	float cv2 = cos(v/2);
	float sv2v = (sin(v/2)/v);
	float sv2vX = sv2v*sX;
	float sv2vY = sv2v*sY;
	float sv2vZ = sv2v*sZ;
	PGMATRIX_DATA(pga44_2_1,1,1) = cv2;
	PGMATRIX_DATA(pga44_2_1,1,2) = -sv2vX;
	PGMATRIX_DATA(pga44_2_1,1,3) = -sv2vY;
	PGMATRIX_DATA(pga44_2_1,1,4) = -sv2vZ;
	PGMATRIX_DATA(pga44_2_1,2,1) = sv2vX;
	PGMATRIX_DATA(pga44_2_1,2,2) = cv2;
	PGMATRIX_DATA(pga44_2_1,2,3) = sv2vZ;
	PGMATRIX_DATA(pga44_2_1,2,4) = -sv2vY;
	PGMATRIX_DATA(pga44_2_1,3,1) = -sv2vY;
	PGMATRIX_DATA(pga44_2_1,3,2) = sv2vZ;
	PGMATRIX_DATA(pga44_2_1,3,3) = cv2;
	PGMATRIX_DATA(pga44_2_1,3,4) = -sv2vY;
	PGMATRIX_DATA(pga44_2_1,4,1) = -sv2vZ;
	PGMATRIX_DATA(pga44_2_1,4,2) = -sv2vY;
	PGMATRIX_DATA(pga44_2_1,4,3) = sv2vX;
	PGMATRIX_DATA(pga44_2_1,4,4) = cv2;

	PGMATRIX_MULTIPLY_COPY (pga44_2_2,pga44_2_1,pPant);
	PGMATRIX_TRANSPOSE_COPY (pga44_2_1,pga44_2_1);
	PGMATRIX_MULTIPLY_COPY (pga44_2_2,pga44_2_2,pga44_2_1);
	PGMATRIX_ADD_COPY (pP,pga44_2_2,pQa);
}

void attitudeCorrection (PGMATRIX pQ, PGMATRIX pP, PGMATRIX pRa, PGMATRIX pQmeas)
{
	float norm = quatNorm(pQ);
	float normInvMin = (-1/norm);
	
	PGMATRIX_ZEROES (pga54_1_1);
	PGMATRIX_DATA(pga54_1_1,1,1) = 1;
	PGMATRIX_DATA(pga54_1_1,2,2) = 1;
	PGMATRIX_DATA(pga54_1_1,3,3) = 1;
	PGMATRIX_DATA(pga54_1_1,4,4) = 1;
	PGMATRIX_DATA(pga54_1_1,5,1) = PGMATRIX_DATA(pQ,1,1)*normInvMin;
	PGMATRIX_DATA(pga54_1_1,5,2) = PGMATRIX_DATA(pQ,2,1)*normInvMin;
	PGMATRIX_DATA(pga54_1_1,5,3) = PGMATRIX_DATA(pQ,3,1)*normInvMin;
	PGMATRIX_DATA(pga54_1_1,5,4) = PGMATRIX_DATA(pQ,4,1)*normInvMin;

	PGMATRIX_DATA(pga51_1_1,1,1) = PGMATRIX_DATA(pQ,1,1);
	PGMATRIX_DATA(pga51_1_1,2,1) = PGMATRIX_DATA(pQ,2,1);
	PGMATRIX_DATA(pga51_1_1,3,1) = PGMATRIX_DATA(pQ,3,1);
	PGMATRIX_DATA(pga51_1_1,4,1) = PGMATRIX_DATA(pQ,4,1);
	PGMATRIX_DATA(pga51_1_1,5,1) = 1 - norm;

	PGMATRIX_DATA(pga51_1_2,1,1) = PGMATRIX_DATA(pQmeas,1,1);
	PGMATRIX_DATA(pga51_1_2,2,1) = PGMATRIX_DATA(pQmeas,2,1);
	PGMATRIX_DATA(pga51_1_2,3,1) = PGMATRIX_DATA(pQmeas,3,1);
	PGMATRIX_DATA(pga51_1_2,4,1) = PGMATRIX_DATA(pQmeas,4,1);
	PGMATRIX_DATA(pga51_1_2,5,1) = 0;

	// K = Ppred*C'*inv(C*Ppred*C' + R)
	PGMATRIX_TRANSPOSE_COPY (pga55_1_2, pga54_1_1);
	PGMATRIX_MULTIPLY_COPY (pga45_1_2, pP, pga55_1_2);
	PGMATRIX_MULTIPLY_COPY (pga55_1_1, pga54_1_1, pga45_1_2);
	PGMATRIX_ADD (pga55_1_1, pRa);
	PGMATRIX_INVERSE (pga55_1_1);
	PGMATRIX_TRIPLEMULTIPLY_COPY (pga45_1_1, pP, pga55_1_2, pga55_1_1, pga45_1_2);

    // Q = Qpred + K*(QmeasAmp - H)
	PGMATRIX_SUBSTRACT (pga51_1_2, pga51_1_1);
	PGMATRIX_MULTIPLY_COPY (pga41_1_1, pga45_1_1, pga51_1_2);
	PGMATRIX_ADD (pQ, pga41_1_1);

	// P = (eye(4) - K*C)*Ppred
	PGMATRIX_MULTIPLY_COPY (pga44_1_1, pga45_1_1, pga54_1_1);
	PGMATRIX_SUBSTRACT_COPY (pga44_1_1, pga44_1_2, pga44_1_1);
	PGMATRIX_MULTIPLY_COPY (pP, pga44_1_1, pP);
}

//-------------------------------------------------
// POSITION MATHS
//-------------------------------------------------

void inertialVelPropagation (PGMATRIX pV, PGMATRIX pQ, PGMATRIX pAcelerometro, PGMATRIX initialGravity, char gravity)
{
	// aqui checamos se usamos gravity como sendo [0 0 g]' ou o valor inicial medido
	if (gravity == 'i')
	{
		PGMATRIX_MULTIPLY_CONST_COPY(pga31_1_2, initialGravity, (-1));
	}
	else
	{
		if (gravity == 'n')
		{
			PGMATRIX_DATA(pga31_1_2,1,1) = 0;
			PGMATRIX_DATA(pga31_1_2,2,1) = 0;
			PGMATRIX_DATA(pga31_1_2,3,1) = GRAVITY;
		}
	}

	// integracao euler normal
	quatTransformation(pga31_1_1,pQ,pAcelerometro);
	PGMATRIX_MULTIPLY_CONST(pga31_1_1,(float)TA_IN_S);
	PGMATRIX_MULTIPLY_CONST(pga31_1_2,(float)TA_IN_S);
	PGMATRIX_ADD_COPY(pga31_1_3,pga31_1_1,pga31_1_2);
	PGMATRIX_ADD(pV,pga31_1_3);
}

void inertialPosPropagation (PGMATRIX pR, PGMATRIX pV)
{
	// integracao euler normal
	PGMATRIX_MULTIPLY_CONST_COPY(pga31_1_1,pV,(float)TA_IN_S);
	PGMATRIX_ADD(pR,pga31_1_1);
}

*/